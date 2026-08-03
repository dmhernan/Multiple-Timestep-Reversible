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

static int pl_contains(PairList *pl, int i, int j)
{
    int k;
    if (i > j) { int tmp = i; i = j; j = tmp; }
    for (k = 0; k < pl->n; k++)
        if (pl->p[k].i == i && pl->p[k].j == j) return 1;
    return 0;
}

static void pl_unique(PairList *pl)
{
    PairList tmp;
    int k;
    pl_init(&tmp, pl->n > 0 ? pl->n : 1);
    for (k = 0; k < pl->n; k++)
        if (!pl_contains(&tmp, pl->p[k].i, pl->p[k].j))
            pl_append(&tmp, pl->p[k].i, pl->p[k].j);
    pl_free(pl);
    *pl = tmp;
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
    int k, b;
    int *keplevels = (int *)calloc(n, sizeof(int));
    if (!keplevels) { fprintf(stderr, "calc_pairs: malloc failed\n"); exit(1); }

    pl_init(indv, arru->n > 0 ? arru->n : 1);
    for (k = 0; k < arru->n; k++) {
        int i = arru->p[k].i;
        int j = arru->p[k].j;
        if (levc[i][j] == lev)
            pl_append(indv, i, j);
    }

    for (k = 0; k < arru->n; k++) {
        int i   = arru->p[k].i;
        int j   = arru->p[k].j;
        int lij = levc[i][j];
        if (lij > keplevels[i]) keplevels[i] = lij;
        if (lij > keplevels[j]) keplevels[j] = lij;
    }

    *nk = 0;
    for (b = 0; b < n; b++)
        if (keplevels[b] == lev)
            kpairs[(*nk)++] = b;

    free(keplevels);
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
    int *id_repeat, n_repeat;
    int is_final;
    (void)n;

    if (indv->n == 0) return 0;

    is_final = (ind == hsub);

    id_repeat = (int *)malloc(2 * arru->n * sizeof(int));
    if (!id_repeat) { fprintf(stderr, "update_levels: malloc failed\n"); exit(1); }

    pl_init(&indv_pos, indv->n);
    pl_init(&indv_neg, indv->n);
    any_pos = 0; any_neg = 0;

    for (k = 0; k < indv->n; k++) {
        int i = indv->p[k].i;
        int j = indv->p[k].j;
        int D = levcp[k] - levc[i][j];
        if (D > 0) { pl_append(&indv_pos, i, j); any_pos = 1; }
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
            for (k = 0; k < indv_pos.n; k++) {
                int i = indv_pos.p[k].i;
                int j = indv_pos.p[k].j;
                int p2;
                for (p2 = 0; p2 < indv->n; p2++) {
                    if (indv->p[p2].i == i && indv->p[p2].j == j) {
                        levc[i][j] = levcp[p2];
                        levc[j][i] = levcp[p2];
                        break;
                    }
                }
            }
            pl_free(&indv_pos); pl_free(&indv_neg);
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
        pl_unique(&arru_new);
        pl_free(arru);
        *arru = arru_new;
        free(inR);

        n_repeat = 0;
        for (k = 0; k < arru->n; k++) {
            int i = arru->p[k].i;
            int j = arru->p[k].j;
            int found_i = 0, found_j = 0, b2;
            for (b2 = 0; b2 < n_repeat; b2++) {
                if (id_repeat[b2] == i) found_i = 1;
                if (id_repeat[b2] == j) found_j = 1;
            }
            if (!found_i) id_repeat[n_repeat++] = i;
            if (!found_j) id_repeat[n_repeat++] = j;
        }

        /* Rare edge case: simultaneous disjoint demotion */
        if (is_final && any_neg) {
            for (k = 0; k < indv_neg.n; k++) {
                int i = indv_neg.p[k].i;
                int j = indv_neg.p[k].j;
                int in_repeat_i = 0, in_repeat_j = 0, b2;
                for (b2 = 0; b2 < n_repeat; b2++) {
                    if (id_repeat[b2] == i) in_repeat_i = 1;
                    if (id_repeat[b2] == j) in_repeat_j = 1;
                }
                if (!in_repeat_i && !in_repeat_j) {
                    levc[i][j]--;
                    levc[j][i]--;
                    G_RARE_EDGE++;
                    if (G_RARE_FILE)
                        fprintf(G_RARE_FILE, "%.6f  %s  %d\n",
                                G_TIME / 365.25,
                                G_RARE_FWD ? "fwd" : "bwd",
                                G_RARE_N);
                }
            }
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
        for (k = 0; k < indv_pos.n; k++) {
            int i = indv_pos.p[k].i;
            int j = indv_pos.p[k].j;
            /* Find index of this pair in indv to get levcp value */
            int p2;
            for (p2 = 0; p2 < indv->n; p2++) {
                if (indv->p[p2].i == i && indv->p[p2].j == j) {
                    levc[i][j] = levcp[p2];
                    levc[j][i] = levcp[p2];
                    break;
                }
            }
        }

        pl_free(&indv_pos); pl_free(&indv_neg);
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
        }
    }

    pl_free(&indv_pos); pl_free(&indv_neg);
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
    double **Q0     = alloc2d(3, n);
    double **v0     = alloc2d(3, n);
    int    *kpairs  = (int *)malloc(n * sizeof(int));
    int nk;

    if (!kpairs) { fprintf(stderr, "driftop: malloc failed\n"); exit(1); }

    /* Refine arru: keep only pairs at lev or above */
    pl_init(&arru, arru_in->n > 0 ? arru_in->n : 1);
    for (k = 0; k < arru_in->n; k++) {
        int pi = arru_in->p[k].i;
        int pj = arru_in->p[k].j;
        if (levc[pi][pj] >= lev)
            pl_append(&arru, pi, pj);
    }
    pl_unique(&arru);

    for (i = 1; i <= hsub; i++) {
        int loop = 0;
        pl_init(&arru_saved, arru.cap);
        pl_copy(&arru_saved, &arru);
        repeated = 0;

        /* Save state */
        for (k = 0; k < 3; k++) {
            int b;
            for (b = 0; b < n; b++) {
                Q0[k][b] = Q[k][b];
                v0[k][b] = v[k][b];
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

            levmax = 0;
            {
                int pi, pj;
                for (pi = 0; pi < n; pi++)
                    for (pj = 0; pj < n; pj++)
                        if (levc[pi][pj] > levmax) levmax = levc[pi][pj];
            }

            if (lev < levmax)
                driftop(Q, v, m, n, lev+1, hlev, levc, rlev, rlev_sq, &arru);

            interact(Q, v, hlev[lev-1] / 2.0, m, n, &indv);

            /* Allocate small 1D levcp for this indv */
            levcp = (int *)malloc(indv.n * sizeof(int));
            if (!levcp) { fprintf(stderr, "driftop: levcp malloc failed\n"); exit(1); }

            calc_levels_pairs(Q, v, m, rlev, rlev_sq, hlev[0], hlev, levc, levcp, &indv);

            flag = update_levels(Q0, v0, Q, v, levc, levcp,
                                 i, hsub, &indv, lev, &arru, n, 0);

            free(levcp);
            pl_free(&indv);

            if (flag) {
                repeated = 1;
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
                    pl_unique(&tmp);
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
    }

    pl_free(&arru);
    free2d(Q0, 3);
    free2d(v0, 3);
    free(kpairs);
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
            int pi, pj;
            for (pi = 0; pi < n; pi++)
                for (pj = 0; pj < n; pj++)
                    if (levc[pi][pj] > levmax) levmax = levc[pi][pj];
        }

        if (levmax > 1)
            driftop(Q, v, m, n, 2, hlev, levc, rlev, rlev_sq, &arru);

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
    free2d(Q0, 3);
    free2d(v0, 3);
    free(kpairs);
}
