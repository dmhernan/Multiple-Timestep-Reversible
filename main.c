#include "mtr.h"
#include <sys/stat.h>

int    G_FLAG    = 0;
double G_RHILL   = 0.0;
int    G_STEP    = 0;
long **G_SUBSTEPS    = NULL;
long **G_REPEATS     = NULL;
long **G_REPEATS_SELF = NULL;
long   G_GLOBAL_STEPS   = 0;
long   G_GLOBAL_REPEATS = 0;
long   G_GLOBAL_FORCED  = 0;
FILE  *G_GLOBAL_REP_FILE = NULL;
int  **G_MAX_LEV     = NULL;
double **G_MAX_LEV_T = NULL;
double **G_MAX_LEV_SEP = NULL;
int    G_N           = 0;
long   G_RARE_EDGE   = 0;
FILE  *G_RARE_FILE   = NULL;
double G_TIME        = 0.0;
int    G_RARE_FWD    = 1;
int    G_RARE_N      = 0;

/* Compute separation between DHC bodies j and k, scaled by G_RHILL */
static double sep_rhill(double **Q, int j, int k)
{
    double d, s = 0.0;
    int c;
    for (c = 0; c < 3; c++) {
        d = Q[c][j] - Q[c][k];
        s += d * d;
    }
    return sqrt(s) / G_RHILL;
}

int main(int argc, char *argv[])
{
    int    n, nstep, steps[LEVMAX];
    double h, tmax, t;
    int    output_interval;
    int    use_testcase;
    double E0, E, L0[3], L[3];
    double a1, em1, a2, em2;
    FILE  *fcons;
    double rlev[LEVMAX], hlev[LEVMAX], rlev_sq[LEVMAX];

    double *m   = NULL;
    double dt_snap = 0.0;
    char   snap_dir[256] = "snap";
    int    snap_every = 0;
    double **q  = NULL, **v  = NULL;
    double **Q  = NULL, **P  = NULL;
    int    **levc = NULL;

    use_testcase = 0;
    if (argc > 1) use_testcase = atoi(argv[1]);

    if (use_testcase == 1) {
        h               = 0.01 * 365.25;
        tmax            = 100.0 * 365.25;
        output_interval = 1;
        G_FLAG          = 0;
        G_RHILL         = 0.0;
        n = 5;
        m    = (double *)malloc(n * sizeof(double));
        q    = alloc2d(3, n); v    = alloc2d(3, n);
        Q    = alloc2d(3, n); P    = alloc2d(3, n);
        levc = alloc2d_int(n, n);
        in3boddunctwo(m, q, v, &n);
        initv(h, 30.0, 2.0, 3, 30, rlev, hlev, steps);

    } else if (use_testcase == 2) {
        h               = 0.03 * 365.25;
        tmax            = 3000.0 * 365.25;
        output_interval = 1;
        G_FLAG          = 1;
        G_RHILL         = 0.436 * pow(50.0, 1.0/3.0);
        n = 5;
        m    = (double *)malloc(n * sizeof(double));
        q    = alloc2d(3, n); v    = alloc2d(3, n);
        Q    = alloc2d(3, n); P    = alloc2d(3, n);
        levc = alloc2d_int(n, n);
        insolardunc(m, q, v, &n);
        initv(h, 3.0, 2.0, 4, 30, rlev, hlev, steps);

    } else if (use_testcase == 3) {
        /* Same as testcase 2 but using Q_min/r_Hill criterion (flag=2) */
        h               = 0.03 * 365.25;
        tmax            = 3000.0 * 365.25;
        output_interval = 1;
        G_FLAG          = 2;
        G_RHILL         = 0.436 * pow(50.0, 1.0/3.0);
        n = 5;
        m    = (double *)malloc(n * sizeof(double));
        q    = alloc2d(3, n); v    = alloc2d(3, n);
        Q    = alloc2d(3, n); P    = alloc2d(3, n);
        levc = alloc2d_int(n, n);
        insolardunc(m, q, v, &n);
        initv(h, 3.0, 2.0, 4, 30, rlev, hlev, steps);

    } else {
        SimState sim;
        int i;
        const char *fname = argc > 1 ? argv[1] : "input.txt";

        /* Pre-scan for n_bodies before allocating */
        {
            FILE *fp = fopen(fname, "r");
            char line[512], key[64];
            double val;
            if (!fp) { fprintf(stderr, "cannot open '%s'\n", fname); return 1; }
            n = 0;
            while (fgets(line, sizeof(line), fp)) {
                if (line[0] == '#' || line[0] == '\n') continue;
                if (sscanf(line, "%s %lf", key, &val) == 2)
                    if (strcmp(key, "n_bodies") == 0) { n = (int)val; break; }
            }
            fclose(fp);
            if (n <= 0) { fprintf(stderr, "n_bodies not found in '%s'\n", fname); return 1; }
        }

        sim.n    = n;
        sim.m    = (double *)malloc(n * sizeof(double));
        sim.Q    = alloc2d(3, n);
        sim.v    = alloc2d(3, n);
        sim.levc = alloc2d_int(n, n);
        q    = alloc2d(3, n);
        v    = alloc2d(3, n);
        read_input(fname, &sim, q, v);
        n               = sim.n;
        h               = sim.h;
        tmax            = sim.tmax;
        output_interval = sim.output_interval;
        G_FLAG          = sim.flag;
        G_RHILL         = sim.rhill;
        dt_snap         = sim.dt_snap;
        strcpy(snap_dir, sim.snap_dir);
        m    = sim.m;
        Q    = alloc2d(3, n); P = alloc2d(3, n);
        levc = alloc2d_int(n, n);
        for (i = 0; i < LEVMAX; i++) { rlev[i] = sim.rlev[i]; hlev[i] = sim.hlev[i]; }
    }

    convert_cart(m, n, q, v, Q, P);

    /* Precompute rlev_sq[i] = rlev[i]^2 for sqrt-free level comparisons */
    {
        int i;
        for (i = 0; i < LEVMAX; i++)
            rlev_sq[i] = rlev[i] * rlev[i];
    }

    calc_levels(Q, v, m, n, rlev, h, levc);
    consqv(m, n, Q, v, &E0, L0);

    /* Allocate per-pair counters */
    G_N           = n;
    G_SUBSTEPS    = alloc2d_long(n, n);
    G_REPEATS     = alloc2d_long(n, n);
    G_REPEATS_SELF = alloc2d_long(n, n);
    G_MAX_LEV     = alloc2d_int(n, n);
    G_MAX_LEV_T   = alloc2d(n, n);
    G_MAX_LEV_SEP = alloc2d(n, n);

    fcons = fopen("data/fcons.txt", "w");
    if (!fcons) fcons = fopen("fcons.txt", "w");
    if (!fcons) { fprintf(stderr, "cannot open fcons.txt\n"); return 1; }

    G_GLOBAL_REP_FILE = fopen("data/global_repeats.txt", "w");
    if (!G_GLOBAL_REP_FILE) G_GLOBAL_REP_FILE = fopen("global_repeats.txt", "w");
    if (G_GLOBAL_REP_FILE)
        fprintf(G_GLOBAL_REP_FILE, "# step  time_yr  repetitions  forced\n");

    /* Initial output */
    if (use_testcase == 1) {
        calcorb(Q, v, m, 1, 2, &a1, &em1);
        calcorb(Q, v, m, 3, 4, &a2, &em2);
        fprintf(fcons, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
                0.0, E0, a1, a2, em1, em2,
                (double)levc[1][2], (double)levc[3][4]);
    } else if (use_testcase == 2 || use_testcase == 3) {
        fprintf(fcons, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %d %d %d %d %d %d\n",
                0.0, E0,
                sep_rhill(Q,1,2), sep_rhill(Q,1,3), sep_rhill(Q,1,4),
                sep_rhill(Q,2,3), sep_rhill(Q,2,4), sep_rhill(Q,3,4),
                levc[1][2], levc[1][3], levc[1][4],
                levc[2][3], levc[2][4], levc[3][4]);
    } else {
        fprintf(fcons, "%.16e %.16e\n", 0.0, E0);
    }

    /* Snapshot setup (general input-file case only) */
    if (dt_snap > 0.0) {
        snap_every = (int)(dt_snap / h + 0.5);
        if (snap_every < 1) snap_every = 1;
        mkdir(snap_dir, 0755);
        write_snapshot(snap_dir, 0, 0.0, m, Q, v, n);
    }

    /* Main integration loop */
    t     = 0.0;
    nstep = 0;

    while (t < tmax) {
        nstep++;
        G_STEP = nstep;
        G_TIME = t;
        mtr_step(Q, v, m, n, levc, hlev, rlev, rlev_sq);
        t += h;

        if (nstep % output_interval == 0) {
            consqv(m, n, Q, v, &E, L);
            if (use_testcase == 1) {
                calcorb(Q, v, m, 1, 2, &a1, &em1);
                calcorb(Q, v, m, 3, 4, &a2, &em2);
                fprintf(fcons, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
                        t, E, a1, a2, em1, em2,
                        (double)levc[1][2], (double)levc[3][4]);
            } else if (use_testcase == 2 || use_testcase == 3) {
                fprintf(fcons, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %d %d %d %d %d %d\n",
                        t, E,
                        sep_rhill(Q,1,2), sep_rhill(Q,1,3), sep_rhill(Q,1,4),
                        sep_rhill(Q,2,3), sep_rhill(Q,2,4), sep_rhill(Q,3,4),
                        levc[1][2], levc[1][3], levc[1][4],
                        levc[2][3], levc[2][4], levc[3][4]);
            } else {
                fprintf(fcons, "%.16e %.16e\n", t, E);
            }
        }

        if (snap_every > 0 && nstep % snap_every == 0)
            write_snapshot(snap_dir, nstep / snap_every, t, m, Q, v, n);
    }

    /* Final diagnostics */
    consqv(m, n, Q, v, &E, L);
    fprintf(stderr, "\nnstep    = %d\n", nstep);
    fprintf(stderr, "dE/E     = %.16e\n", (E   - E0)    / E0);
    fprintf(stderr, "dLz/Lz   = %.16e\n", (L[2] - L0[2]) / L0[2]);

    /* Substep counter summary */
    {
        long total_substeps = 0, total_repeats = 0;
        int i, j;
        for (i = 1; i < n; i++)
            for (j = i+1; j < n; j++) {
                total_substeps += G_SUBSTEPS[i][j];
                total_repeats  += G_REPEATS[i][j];
            }
        fprintf(stderr, "substeps        = %ld\n", total_substeps);
        fprintf(stderr, "repeats         = %ld\n", total_repeats);
        fprintf(stderr, "repeat_fraction = %.6f\n",
                total_substeps > 0 ? (double)total_repeats / total_substeps : 0.0);
        fprintf(stderr, "rare_edge       = %ld\n", G_RARE_EDGE);
        fprintf(stderr, "global_steps    = %ld\n", G_GLOBAL_STEPS);
        fprintf(stderr, "global_repeats  = %ld\n", G_GLOBAL_REPEATS);
        fprintf(stderr, "global_repeat_fraction = %.6f\n",
                G_GLOBAL_STEPS > 0 ? (double)G_GLOBAL_REPEATS / G_GLOBAL_STEPS : 0.0);
        fprintf(stderr, "global_forced   = %ld  (global steps accepted at the %d-repetition cap)\n",
                G_GLOBAL_FORCED, MAX_GLOBAL_REPEAT);
        fprintf(stderr, "global_forced_fraction = %.6f\n",
                G_GLOBAL_STEPS > 0 ? (double)G_GLOBAL_FORCED / G_GLOBAL_STEPS : 0.0);
    }

    /* Max substep level per pair (testcase 2 only) */
    if (use_testcase == 2 || use_testcase == 3) {
        const char *names[5] = {"Jupiter","Saturn","Uranus","Neptune"};
        int i, j;
        fprintf(stderr, "\nMax substep level per pair:\n");
        for (i = 1; i < n; i++)
            for (j = i+1; j < n; j++)
                fprintf(stderr, "  %s--%s: level %d at t = %.2f yr, sep = %.4f au\n",
                        names[i-1], names[j-1],
                        G_MAX_LEV[i][j],
                        G_MAX_LEV_T[i][j] / 365.25,
                        G_MAX_LEV_SEP[i][j]);
    }

    /* Write per-pair counter file */
    {
        FILE *fcnt = fopen("data/substep_counts.txt", "w");
        int i, j;
        if (!fcnt) fcnt = fopen("substep_counts.txt", "w");
        if (fcnt) {
            fprintf(fcnt, "# i j substeps repeats repeats_self fraction\n");
            for (i = 1; i < n; i++)
                for (j = i+1; j < n; j++) {
                    long s  = G_SUBSTEPS[i][j];
                    long r  = G_REPEATS[i][j];
                    long rs = G_REPEATS_SELF[i][j];
                    fprintf(fcnt, "%d %d %ld %ld %ld %.6f\n",
                            i, j, s, r, rs,
                            s > 0 ? (double)r / s : 0.0);
                }
            fclose(fcnt);
        }
    }

    fclose(fcons);
    if (G_GLOBAL_REP_FILE) fclose(G_GLOBAL_REP_FILE);

    free(m);
    free2d(q, 3); free2d(v, 3);
    free2d(Q, 3); free2d(P, 3);
    free2d_int(levc, n);
    free2d_long(G_SUBSTEPS, n);
    free2d_long(G_REPEATS, n);
    free2d_long(G_REPEATS_SELF, n);
    free2d_int(G_MAX_LEV, n);
    free2d(G_MAX_LEV_T, n);
    free2d(G_MAX_LEV_SEP, n);

    return 0;
}
