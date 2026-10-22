#include "mtr.h"

/*
   integrator.c
   Core integration routines for the time-reversible multi-timestep
   integrator.  Based on Hernandez & Dehnen (2024).

   Mixed coordinates throughout:
     Q[3][n]  -- positions in DHC coordinates
     v[3][n]  -- velocities in Cartesian

   Key functions:
     calc_pairs     -- find interaction pairs and Kepler bodies at level
     update_levels  -- handle level promotion/demotion with reversibility
     driftop        -- recursive multi-level drift operator
     mtr_step       -- top-level time-reversible integrator step
*/

/* ───────────────────────────────────────────────────────────────────────────
   Helper: PairList operations (beyond those in mtr.h)
   ─────────────────────────────────────────────────────────────────────────── */

/* ───────────────────────────────────────────────────────────────────────────
   Persistent scratch buffers.
   These replace per-call/per-pass allocations of n-sized arrays.  They are
   (re)allocated only when n changes.  calc_pairs and the body-list helper
   are never re-entered, so one shared set suffices; driftop's saved state
   is held while deeper levels run, so it gets one buffer per level.
   ─────────────────────────────────────────────────────────────────────────── */
/* Pairs demoted since the enclosing level last looked (DEMOTION FIX).
   A pair demoted at the end of a level-(l+1) block is, at that same instant,
   re-examined by the enclosing level-l substep, before it receives any further
   kick; if that substep also ends its block, it may be demoted again. */
static PairList JD = {NULL, 0, 0};
static void jd_push(int i, int j) {
    if (JD.p == NULL) pl_init(&JD, 64);
    pl_append(&JD, i, j);
}
static int     sc_n = -1;
static int    *sc_keplev = NULL;   /* zero between uses                     */
static char   *sc_mark   = NULL;   /* zero between uses                     */
static int    *sc_list   = NULL;   /* bodies touched in the current call    */
static double **sc_Q0[LEVMAX + 2];
static double **sc_v0[LEVMAX + 2];
static int    *sc_kp[LEVMAX + 2];

static void scratch_ensure(int n)
{
    int l;
    if (n == sc_n) return;
    free(sc_keplev); free(sc_mark); free(sc_list);
    for (l = 0; l < LEVMAX + 2; l++) {
        if (sc_Q0[l]) free2d(sc_Q0[l], 3);
        if (sc_v0[l]) free2d(sc_v0[l], 3);
        free(sc_kp[l]);
        sc_Q0[l] = NULL; sc_v0[l] = NULL; sc_kp[l] = NULL;
    }
    sc_keplev = (int *)calloc(n, sizeof(int));
    sc_mark   = (char *)calloc(n, 1);
    sc_list   = (int *)malloc(n * sizeof(int));
    if (!sc_keplev || !sc_mark || !sc_list) {
        fprintf(stderr, "scratch_ensure: malloc failed\n"); exit(1);
    }
    sc_n = n;
}

static void scratch_level(int n, int lev)
{
    scratch_ensure(n);
    if (lev < 0 || lev >= LEVMAX + 2) {
        fprintf(stderr, "scratch_level: level %d out of range\n", lev); exit(1);
    }
    if (!sc_Q0[lev]) {
        sc_Q0[lev] = alloc2d(3, n);
        sc_v0[lev] = alloc2d(3, n);
        sc_kp[lev] = (int *)malloc(n * sizeof(int));
        if (!sc_kp[lev]) { fprintf(stderr, "scratch_level: malloc failed\n"); exit(1); }
    }
}

/* EXPERIMENTAL FIX: after a rewind, pairs that were reset and sat below the
   current level (lev_saved > lev) get their substep-start level back, so that
   level changes made inside the discarded trial do not leak into the repeat.
   Promoted pairs at the current level (lev_saved == lev) keep their new level. */
static void restore_deeper_levels(int **levc, const PairList *saved, const int *lev_saved,
                                  const PairList *rewound, int lev, int n)
{
    int p;
    scratch_ensure(n);
    for (p = 0; p < rewound->n; p++) { sc_mark[rewound->p[p].i] = 1; sc_mark[rewound->p[p].j] = 1; }
    for (p = 0; p < saved->n; p++) {
        int a = saved->p[p].i, b = saved->p[p].j;
        if ((sc_mark[a] || sc_mark[b]) && lev_saved[p] > lev) {
            levc[a][b] = lev_saved[p];
            levc[b][a] = lev_saved[p];
        }
    }
    for (p = 0; p < rewound->n; p++) { sc_mark[rewound->p[p].i] = 0; sc_mark[rewound->p[p].j] = 0; }
}

static void pl_all_pairs(PairList *pl, int n)
{
    int i, j;
    int cap = (n*(n-1))/2;
    if (cap < 1) cap = 1;
    pl_init(pl, cap);
    for (i = 1; i < n; i++)
        for (j = i+1; j < n; j++)
            pl_append(pl, i, j);
}

/* ───────────────────────────────────────────────────────────────────────────
   calc_pairs
   Find interaction pairs (indv) and Kepler bodies (kpairs) at level lev.
   ─────────────────────────────────────────────────────────────────────────── */
void calc_pairs(int lev, int **levc, int n,
                PairList *arru, PairList *indv, int *kpairs, int *nk)
{
    /* One pass over arru builds both indv and the per-body Kepler level.
       Only bodies appearing in arru are touched (previously all n bodies
       were allocated, zeroed and scanned on every pass).  kpairs is filled in
       order of first appearance in arru rather than ascending body index.
       The order does not affect results: each body's Kepler step reads and
       writes only that body's own Q and v, and the solver keeps no state
       between calls. */
    int k, t, nt = 0;
    int *keplevels;
    scratch_ensure(n);
    keplevels = sc_keplev;

    pl_init(indv, arru->n > 0 ? arru->n : 1);
    for (k = 0; k < arru->n; k++) {
        int i   = arru->p[k].i;
        int j   = arru->p[k].j;
        int lij = levc[i][j];
        if (lij == lev)
            pl_append(indv, i, j);
        if (!sc_mark[i]) { sc_mark[i] = 1; sc_list[nt++] = i; }
        if (!sc_mark[j]) { sc_mark[j] = 1; sc_list[nt++] = j; }
        if (lij > keplevels[i]) keplevels[i] = lij;
        if (lij > keplevels[j]) keplevels[j] = lij;
    }

    *nk = 0;
    for (t = 0; t < nt; t++) {
        int b = sc_list[t];
        if (keplevels[b] == lev)
            kpairs[(*nk)++] = b;
        keplevels[b] = 0;      /* restore scratch to all-zero */
        sc_mark[b]   = 0;
    }
}

/* ───────────────────────────────────────────────────────────────────────────
   update_levels
   Handle level promotion and demotion after a substep.
   ─────────────────────────────────────────────────────────────────────────── */
int update_levels(double **Q0, double **v0, double **Q, double **v,
                  int **levc, int *levcp,
                  int ind, int hsub, PairList *indv, int lev,
                  PairList *arru, int n, int force_accept)
{
    int k, b;
    int any_pos, any_neg;
    PairList indv_pos, indv_neg;
    int *pos_idx;   /* pos_idx[q] = index in indv of promoted pair q */
    int *id_repeat, n_repeat;
    int is_final;
    (void)n;

    if (indv->n == 0) return 0;

    is_final = (ind == hsub);

    id_repeat = (int *)malloc(2 * arru->n * sizeof(int));
    if (!id_repeat) { fprintf(stderr, "update_levels: malloc failed\n"); exit(1); }

    pl_init(&indv_pos, indv->n);
    pl_init(&indv_neg, indv->n);
    pos_idx = (int *)malloc(indv->n * sizeof(int));
    if (!pos_idx) { fprintf(stderr, "update_levels: malloc failed\n"); exit(1); }
    any_pos = 0; any_neg = 0;

    for (k = 0; k < indv->n; k++) {
        int i = indv->p[k].i;
        int j = indv->p[k].j;
        int D = levcp[k] - levc[i][j];
        if (D > 0) { pos_idx[indv_pos.n] = k; pl_append(&indv_pos, i, j); any_pos = 1; }
        if (D < 0) { pl_append(&indv_neg, i, j); any_neg = 1; }
    }

    /* Overlapping promotion and demotion: no special action needed.
       The transitive closure below will pull any body shared between
       a promoted and demoted pair into the rewind set, so the demoted
       pair gets repeated and its demotion signal from the failed trial
       is discarded.  It will be re-evaluated after the repeat. */

    /* Handle promotion */
    if (any_pos) {
        if (force_accept) {
            /* Global repetition cap reached: accept the current trial
               state (Qf,vf) exactly as computed, with NO rewind. Still
               update levc bookkeeping for the promoted pairs so later
               steps use the correct (finer) level, but do not request
               a repeat. Signalled to the caller via return value 2,
               distinct from the normal "no repeat needed" case (0). */
            /* pos_idx[k] is this pair's index in indv, recorded when
               indv_pos was built (replaces an O(|indv|) linear search) */
            for (k = 0; k < indv_pos.n; k++) {
                int i = indv_pos.p[k].i;
                int j = indv_pos.p[k].j;
                levc[i][j] = levcp[pos_idx[k]];
                levc[j][i] = levcp[pos_idx[k]];
            }
            pl_free(&indv_pos); pl_free(&indv_neg);
            free(pos_idx);
            free(id_repeat);
            return 2;
        }

        /* Transitive closure: expand rewind set until no new bodies are added.
           inR[b] = 1 if body b is in the rewind set.
           Seed with bodies of promoted pairs, then expand to any body
           connected via a pair with levc >= lev in arru. */
        int *inR = (int *)calloc(n, sizeof(int));
        int changed;
        PairList arru_new;
        if (!inR) { fprintf(stderr, "update_levels: malloc failed\n"); exit(1); }

        /* Global-level shortcut.  At level 1, arru is the complete list of
           all pairs of the n-1 planets.  If every one of them has
           levc >= 1, the pairs form a complete graph, so the closure below
           necessarily contains every planet and arru_new would equal arru,
           pair for pair and in the same order.  In that case skip the
           closure and the copy and keep arru unchanged.  If any pair has
           level 0 (possible only for separations below the smallest cutoff),
           fall through to the general computation. */
        int full = 0;
        if (lev == 1 && n - 1 >= 2 && arru->n == (n - 1) * (n - 2) / 2) {
            full = 1;
            for (k = 0; k < arru->n; k++)
                if (levc[arru->p[k].i][arru->p[k].j] < 1) { full = 0; break; }
        }

        if (!full) {
        /* Seed: bodies of promoted pairs */
        for (k = 0; k < indv_pos.n; k++) {
            inR[indv_pos.p[k].i] = 1;
            inR[indv_pos.p[k].j] = 1;
        }

        /* Expand to fixed point */
        do {
            changed = 0;
            for (k = 0; k < arru->n; k++) {
                int i = arru->p[k].i;
                int j = arru->p[k].j;
                if (levc[i][j] >= lev && (inR[i] || inR[j])) {
                    if (!inR[i]) { inR[i] = 1; changed = 1; }
                    if (!inR[j]) { inR[j] = 1; changed = 1; }
                }
            }
        } while (changed);

        /* Build arru_new: pairs in arru with levc >= lev and at least one
           endpoint in the rewind set */
        pl_init(&arru_new, arru->n > 0 ? arru->n : 1);
        for (k = 0; k < arru->n; k++) {
            int i   = arru->p[k].i;
            int j   = arru->p[k].j;
            int lij = levc[i][j];
            if (lij >= lev && (inR[i] || inR[j]))
                pl_append(&arru_new, i, j);
        }
        /* pl_unique removed: arru_new is a filter of a duplicate-free list */
        pl_free(arru);
        *arru = arru_new;
        } /* end if (!full) */
        free(inR);

        n_repeat = 0;
        {
            char *seen = (char *)calloc(n, 1);
            if (!seen) { fprintf(stderr, "update_levels: malloc failed\n"); exit(1); }
            for (k = 0; k < arru->n; k++) {
                int i = arru->p[k].i;
                int j = arru->p[k].j;
                if (!seen[i]) { seen[i] = 1; id_repeat[n_repeat++] = i; }
                if (!seen[j]) { seen[j] = 1; id_repeat[n_repeat++] = j; }
            }
            /* membership test for the rare-edge block below */
            if (is_final && any_neg) {
                for (k = 0; k < indv_neg.n; k++) {
                    int i = indv_neg.p[k].i;
                    int j = indv_neg.p[k].j;
                    if (!seen[i] && !seen[j]) {
                        levc[i][j]--;
                        levc[j][i]--;
                        jd_push(i, j);
                        G_RARE_EDGE++;
                        if (G_RARE_FILE)
                            fprintf(G_RARE_FILE, "%.6f  %s  %d\n",
                                    G_TIME / 365.25,
                                    G_RARE_FWD ? "fwd" : "bwd",
                                    G_RARE_N);
                    }
                }
            }
            free(seen);
        }

        /* Reset positions and velocities */
        for (k = 0; k < n_repeat; k++) {
            int coord;
            b = id_repeat[k];
            for (coord = 0; coord < 3; coord++) {
                Q[coord][b] = Q0[coord][b];
                v[coord][b] = v0[coord][b];
            }
        }

        /* Update levc for promoted pairs */
        /* pos_idx[k] is this pair's index in indv, recorded when
           indv_pos was built (replaces an O(|indv|) linear search) */
        for (k = 0; k < indv_pos.n; k++) {
            int i = indv_pos.p[k].i;
            int j = indv_pos.p[k].j;
            levc[i][j] = levcp[pos_idx[k]];
            levc[j][i] = levcp[pos_idx[k]];
        }

        pl_free(&indv_pos); pl_free(&indv_neg);
        free(pos_idx);
        free(id_repeat);
        /* Log levc state */
        return 1;
    }

    /* Handle demotion only */
    if (is_final && any_neg) {
        for (k = 0; k < indv_neg.n; k++) {
            int i = indv_neg.p[k].i;
            int j = indv_neg.p[k].j;
            levc[i][j]--;
            levc[j][i]--;
            jd_push(i, j);
        }
    }

    pl_free(&indv_pos); pl_free(&indv_neg);
    free(pos_idx);
    free(id_repeat);
    /* Log levc state */
    return 0;
}

/* ───────────────────────────────────────────────────────────────────────────
   driftop
   Recursive multi-level drift operator.
   ─────────────────────────────────────────────────────────────────────────── */
void driftop(double **Q, double **v, double *m,
             int n, int lev, double *hlev,
             int **levc, double *rlev, double *rlev_sq, PairList *arru_in)
{
    int hsub = (int)round(hlev[0] / hlev[1]);
    int levmax, i, k;
    PairList arru, arru_saved, indv;
    int flag, repeated;
    double **Q0, **v0;
    int    *kpairs;
    int nk;

    /* Per-level persistent buffers (previously allocated on every call) */
    scratch_level(n, lev);
    Q0     = sc_Q0[lev];
    v0     = sc_v0[lev];
    kpairs = sc_kp[lev];

    /* Refine arru: keep only pairs at lev or above */
    pl_init(&arru, arru_in->n > 0 ? arru_in->n : 1);
    for (k = 0; k < arru_in->n; k++) {
        int pi = arru_in->p[k].i;
        int pj = arru_in->p[k].j;
        if (levc[pi][pj] >= lev)
            pl_append(&arru, pi, pj);
    }
    /* pl_unique removed */

    for (i = 1; i <= hsub; i++) {
        int loop = 0;
        pl_init(&arru_saved, arru.cap);
        pl_copy(&arru_saved, &arru);
        repeated = 0;
        int *lev_saved = (int *)malloc((arru.n > 0 ? arru.n : 1) * sizeof(int));
        { int p; for (p = 0; p < arru.n; p++) lev_saved[p] = levc[arru.p[p].i][arru.p[p].j]; }

        /* Save state of the bodies in arru only.  Within this substep only
           those bodies are kicked, drifted (here or at deeper levels), or
           reset by update_levels (its rewind set is a subset of arru), so
           the saved values of all other bodies are never read.  Previously
           all n bodies were copied on every substep. */
        {
            int p, nb = 0, t, b;
            for (p = 0; p < arru.n; p++) {
                int pi = arru.p[p].i, pj = arru.p[p].j;
                if (!sc_mark[pi]) { sc_mark[pi] = 1; sc_list[nb++] = pi; }
                if (!sc_mark[pj]) { sc_mark[pj] = 1; sc_list[nb++] = pj; }
            }
            for (t = 0; t < nb; t++) {
                b = sc_list[t];
                for (k = 0; k < 3; k++) {
                    Q0[k][b] = Q[k][b];
                    v0[k][b] = v[k][b];
                }
                sc_mark[b] = 0;
            }
        }

        while (1) {
            int *levcp;
            loop++;

            calc_pairs(lev, levc, n, &arru, &indv, kpairs, &nk);

            interact(Q, v, hlev[lev-1] / 2.0, m, n, &indv);
            kep_map(Q, v, hlev[lev-1], m, n, kpairs, nk);

            /* Count substeps per pair */
            if (G_SUBSTEPS) {
                int p;
                for (p = 0; p < indv.n; p++) {
                    int pi = indv.p[p].i, pj = indv.p[p].j;
                    if (pi > pj) { int tmp=pi; pi=pj; pj=tmp; }
                    G_SUBSTEPS[pi][pj]++;
                    if (G_MAX_LEV && levc[pi][pj] > G_MAX_LEV[pi][pj]) {
                        double dsep = 0.0; int c;
                        for (c = 0; c < 3; c++) {
                            double dq = Q[c][pi] - Q[c][pj];
                            dsep += dq * dq;
                        }
                        G_MAX_LEV[pi][pj]    = levc[pi][pj];
                        G_MAX_LEV_T[pi][pj]  = G_TIME;
                        G_MAX_LEV_SEP[pi][pj] = sqrt(dsep);
                    }
                }
            }

            /* Deepest level among the ACTIVE pairs only: deeper levels can
               only ever act on pairs in arru, so scanning all n x n entries
               is unnecessary (and costs O(N^2) per substep). */
            levmax = 0;
            {
                int p;
                for (p = 0; p < arru.n; p++) {
                    int l = levc[arru.p[p].i][arru.p[p].j];
                    if (l > levmax) levmax = l;
                }
            }

            int jd_mark = (JD.p != NULL) ? JD.n : 0;
            if (lev < levmax)
                driftop(Q, v, m, n, lev+1, hlev, levc, rlev, rlev_sq, &arru);

            interact(Q, v, hlev[lev-1] / 2.0, m, n, &indv);

            /* DEMOTION FIX: pairs demoted into this level by the deeper call
               (necessarily at the end of this substep) join this substep's
               level check, demotion only.  They are appended after the closing
               kick, so they receive no kick at this level. */
            int n_kicked = indv.n;
            if (JD.p != NULL && JD.n > jd_mark) {
                int q;
                for (q = jd_mark; q < JD.n; q++) {
                    int a = JD.p[q].i, b = JD.p[q].j, p2, dup = 0;
                    if (levc[a][b] != lev) continue;              /* not (or no longer) at this level */
                    for (p2 = n_kicked; p2 < indv.n; p2++)        /* avoid duplicates */
                        if (indv.p[p2].i == a && indv.p[p2].j == b) { dup = 1; break; }
                    if (!dup) pl_append(&indv, a, b);
                }
                JD.n = jd_mark;                                   /* consumed */
            }

            /* Allocate small 1D levcp for this indv */
            levcp = (int *)malloc((indv.n > 0 ? indv.n : 1) * sizeof(int));
            if (!levcp) { fprintf(stderr, "driftop: levcp malloc failed\n"); exit(1); }

            calc_levels_pairs(Q, v, m, rlev, rlev_sq, hlev[0], hlev, levc, levcp, &indv);
            {   /* re-examined pairs may only be demoted, never promoted */
                int q;
                for (q = n_kicked; q < indv.n; q++) if (levcp[q] > lev) levcp[q] = lev;
            }

            flag = update_levels(Q0, v0, Q, v, levc, levcp,
                                 i, hsub, &indv, lev, &arru, n, 0);

            free(levcp);
            pl_free(&indv);

            if (flag) {
                repeated = 1;
                restore_deeper_levels(levc, &arru_saved, lev_saved, &arru, lev, n);
                if (G_REPEATS) {
                    int p;
                    for (p = 0; p < arru.n; p++) {
                        int pi = arru.p[p].i, pj = arru.p[p].j;
                        if (pi > pj) { int tmp=pi; pi=pj; pj=tmp; }
                        G_REPEATS[pi][pj]++;
                        if (G_REPEATS_SELF && levc[pi][pj] > lev)
                            G_REPEATS_SELF[pi][pj]++;
                    }
                }
            } else {
                if (repeated) {
                    PairList tmp;
                    pl_free(&arru);
                    pl_init(&arru, arru_saved.n > 0 ? arru_saved.n : 1);
                    pl_copy(&arru, &arru_saved);
                    /* Remove pairs now below lev */
                    pl_init(&tmp, arru.n > 0 ? arru.n : 1);
                    for (k = 0; k < arru.n; k++) {
                        int pi = arru.p[k].i;
                        int pj = arru.p[k].j;
                        if (levc[pi][pj] >= lev)
                            pl_append(&tmp, pi, pj);
                    }
                    /* pl_unique removed */
                    pl_free(&arru);
                    arru = tmp;
                }
                break;
            }

            if (loop > 50) {
                fprintf(stderr,
                    "driftop: too many repetitions at level %d substep %d\n",
                    lev, i);
                exit(-1);
            }
        }
        pl_free(&arru_saved);
        free(lev_saved);
    }

    pl_free(&arru);
    /* Q0, v0, kpairs are persistent per-level buffers: not freed here */
}

/* ───────────────────────────────────────────────────────────────────────────
   mtr_step
   Top-level time-reversible integrator step.
   ─────────────────────────────────────────────────────────────────────────── */
void mtr_step(double **Q, double **v, double *m,
              int n, int **levc, double *hlev, double *rlev, double *rlev_sq)
{
    double p0[3];
    PairList arru, arru_saved, indv;
    int *kpairs = (int *)malloc(n * sizeof(int));
    int nk;
    int flag, repeated, loop, levmax, k, b;
    int hsub = 2;
    double **Q0 = alloc2d(3, n);
    double **v0 = alloc2d(3, n);

    if (!kpairs) { fprintf(stderr, "mtr_step: malloc failed\n"); exit(1); }

    calc_m(v, m, n, p0);
    map_sun(Q, v, hlev[0] / 2.0, m, n);

    pl_all_pairs(&arru, n);
    pl_init(&arru_saved, arru.cap);
    pl_copy(&arru_saved, &arru);
    repeated = 0;
    int *lev_saved = (int *)malloc((arru.n > 0 ? arru.n : 1) * sizeof(int));
    { int p; for (p = 0; p < arru.n; p++) lev_saved[p] = levc[arru.p[p].i][arru.p[p].j]; }

    for (k = 0; k < 3; k++)
        for (b = 0; b < n; b++) {
            Q0[k][b] = Q[k][b];
            v0[k][b] = v[k][b];
        }

    loop = 0;
    G_GLOBAL_STEPS++;
    {
    int forced = 0;
    while (1) {
        int *levcp;
        int force_this;
        loop++;
        force_this = (loop >= MAX_GLOBAL_REPEAT);

        calc_pairs(1, levc, n, &arru, &indv, kpairs, &nk);

        interact(Q, v, hlev[0] / 2.0, m, n, &indv);
        kep_map(Q, v, hlev[0], m, n, kpairs, nk);

        /* Count substeps per pair */
        if (G_SUBSTEPS) {
            int p;
            for (p = 0; p < indv.n; p++) {
                int pi = indv.p[p].i, pj = indv.p[p].j;
                if (pi > pj) { int tmp=pi; pi=pj; pj=tmp; }
                G_SUBSTEPS[pi][pj]++;
                if (G_MAX_LEV && levc[pi][pj] > G_MAX_LEV[pi][pj]) {
                    double dsep = 0.0; int c;
                    for (c = 0; c < 3; c++) {
                        double dq = Q[c][pi] - Q[c][pj];
                        dsep += dq * dq;
                    }
                    G_MAX_LEV[pi][pj]     = levc[pi][pj];
                    G_MAX_LEV_T[pi][pj]   = G_TIME;
                    G_MAX_LEV_SEP[pi][pj] = sqrt(dsep);
                }
            }
        }

        levmax = 0;
        {
            int p;
            for (p = 0; p < arru.n; p++) {
                int l = levc[arru.p[p].i][arru.p[p].j];
                if (l > levmax) levmax = l;
            }
        }

        if (levmax > 1)
            driftop(Q, v, m, n, 2, hlev, levc, rlev, rlev_sq, &arru);
        if (JD.p != NULL) JD.n = 0;   /* DEMOTION FIX: nothing above the global level */

        interact(Q, v, hlev[0] / 2.0, m, n, &indv);

        /* Allocate small 1D levcp for this indv */
        levcp = (int *)malloc(indv.n * sizeof(int));
        if (!levcp) { fprintf(stderr, "mtr_step: levcp malloc failed\n"); exit(1); }

        calc_levels_pairs(Q, v, m, rlev, rlev_sq, hlev[0], hlev, levc, levcp, &indv);

        flag = update_levels(Q0, v0, Q, v, levc, levcp,
                             0, hsub, &indv, 1, &arru, n, force_this);

        free(levcp);
        pl_free(&indv);

        if (flag == 1) {
            repeated = 1;
            restore_deeper_levels(levc, &arru_saved, lev_saved, &arru, 1, n);
            G_GLOBAL_REPEATS++;
            if (G_REPEATS) {
                int p;
                for (p = 0; p < arru.n; p++) {
                    int pi = arru.p[p].i, pj = arru.p[p].j;
                    if (pi > pj) { int tmp=pi; pi=pj; pj=tmp; }
                    G_REPEATS[pi][pj]++;
                    if (G_REPEATS_SELF && levc[pi][pj] > 1)
                        G_REPEATS_SELF[pi][pj]++;
                }
            }
        } else {
            if (flag == 2) {
                /* Cap reached: trial state (Qf,vf) was accepted as-is,
                   with no further rewind/repeat. */
                forced = 1;
                G_GLOBAL_FORCED++;
            }
            if (repeated) {
                pl_free(&arru);
                pl_init(&arru, arru_saved.cap);
                pl_copy(&arru, &arru_saved);
            }
            break;
        }

        /* Defensive backstop only -- force_this at MAX_GLOBAL_REPEAT
           guarantees the loop always terminates by then, so this should
           never actually trigger. Kept as a safety net against a future
           logic error rather than an expected code path. */
        if (loop > 1000) {
            fprintf(stderr, "mtr_step: unexpected runaway repetition loop (>1000)\n");
            exit(-1);
        }
    }

    if (G_GLOBAL_REP_FILE)
        fprintf(G_GLOBAL_REP_FILE, "%d  %.6f  %d  %d\n",
                G_STEP, G_TIME / 365.25, loop - 1, forced);
    }

    map_sun(Q, v, hlev[0] / 2.0, m, n);
    adjust_sun(v, m, n, p0);

    pl_free(&arru);
    pl_free(&arru_saved);
    free(lev_saved);
    free2d(Q0, 3);
    free2d(v0, 3);
    free(kpairs);
}
