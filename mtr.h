#ifndef MTR_H
#define MTR_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ─── Physical constants ─────────────────────────────────────────────────── */
#define GNEWT  2.95912208232212840e-4  /* G in au^3 day^-2 solar_mass^-1     */

/* ─── Global simulation flags (set once at startup in main.c) ────────────── */
extern int    G_FLAG;    /* 0 = timescale criterion, 1 = Hill radius          */
extern double G_RHILL;   /* Hill radius (au), used if G_FLAG == 1             */
extern int    G_STEP;    /* current step number, for diagnostic output        */
extern long **G_SUBSTEPS;
extern long **G_REPEATS;
extern long **G_REPEATS_SELF;
extern long   G_GLOBAL_STEPS;    /* number of global steps attempted        */
extern long   G_GLOBAL_REPEATS;  /* number of global steps that were repeated */
extern long   G_GLOBAL_FORCED;   /* number of global steps force-accepted at
                                     the MAX_GLOBAL_REPEAT cap              */
extern FILE  *G_GLOBAL_REP_FILE; /* per-global-step repetition count log, or
                                     NULL to disable (set by caller)       */
extern int  **G_MAX_LEV;      /* max substep level reached per pair */
extern double **G_MAX_LEV_T;  /* time (days) when max substep level was reached */
extern double **G_MAX_LEV_SEP; /* separation (au) when max substep level was reached */
extern int    G_N;
extern long   G_RARE_EDGE;

/* Rare event logging */
extern FILE  *G_RARE_FILE;   /* output file, set by caller            */
extern double G_TIME;        /* current physical time (days)          */
extern int    G_RARE_FWD;    /* 1 = forward pass, 0 = backward pass   */
extern int    G_RARE_N;      /* current n value being tested          */

/* ─── Solver status codes ────────────────────────────────────────────────── */
#define SUCCESS  1
#define FAILURE  0

/* ─── Array size limits ──────────────────────────────────────────────────── */
#define NMAX     5000   /* maximum number of bodies (sanity check only)       */
#define LEVMAX   30     /* maximum number of timestep levels                  */
#define PAIRMAX  (500*499/2)  /* max pairs for sanity; resized at runtime     */

/* Cap on global-level repetitions (mtr_step's retry loop) before the trial
   step is forced to be accepted as-is, rather than rewound and repeated
   again. Applies ONLY to the global level; driftop's finer-level substep
   loop keeps its own separate (much larger) hard cap. */
#define MAX_GLOBAL_REPEAT 10

/* ─── State vector for one body ──────────────────────────────────────────── */
typedef struct {
    double x, y, z;      /* position  (au)                */
    double xd, yd, zd;   /* velocity  (au/day)            */
} State;

/* ─── Global simulation state ────────────────────────────────────────────── */
typedef struct {
    int    n;                    /* number of bodies                          */
    double *m;                   /* masses (solar masses)         [n]         */
    double **Q;                  /* DHC positions                 [3][n]      */
    double **v;                  /* velocities                    [3][n]      */
    int    **levc;               /* level assignments             [n][n]      */
    double hlev[LEVMAX];         /* timestep at each level                    */
    double rlev[LEVMAX];         /* radius threshold at each level            */
    int    hsub;                 /* substep multiplier                        */
    int    rsub;                 /* radius ratio between levels               */
    double h;                    /* global timestep (days)                    */
    double t;                    /* current time (days)                       */
    double tmax;                 /* maximum time (days)                       */
    int    output_interval;      /* steps between output writes               */
    int    flag;                 /* 0 = timescale criterion, 1 = Hill radius  */
    double rhill;                /* Hill radius (au), used if flag==1         */
    int    cas;                  /* coordinate case (1 = Dem. Heliocentric)   */
    double dt_snap;              /* snapshot interval; 0 = no snapshots       */
    char   snap_dir[256];        /* snapshot output directory                 */
} SimState;

/* ─── Pair list ──────────────────────────────────────────────────────────── */
typedef struct {
    int i, j;
} Pair;

typedef struct {
    Pair  *p;   /* heap-allocated array of pairs                              */
    int    n;   /* number of pairs currently stored                           */
    int    cap; /* allocated capacity                                         */
} PairList;

/* ─── Heap allocation helpers ────────────────────────────────────────────── */
static inline double **alloc2d(int rows, int cols) {
    int i;
    double **a = (double **)malloc(rows * sizeof(double *));
    if (!a) { fprintf(stderr, "alloc2d: out of memory\n"); exit(1); }
    for (i = 0; i < rows; i++) {
        a[i] = (double *)calloc(cols, sizeof(double));
        if (!a[i]) { fprintf(stderr, "alloc2d: out of memory\n"); exit(1); }
    }
    return a;
}

static inline void free2d(double **a, int rows) {
    int i;
    for (i = 0; i < rows; i++) free(a[i]);
    free(a);
}

static inline int **alloc2d_int(int rows, int cols) {
    int i;
    int **a = (int **)malloc(rows * sizeof(int *));
    if (!a) { fprintf(stderr, "alloc2d_int: out of memory\n"); exit(1); }
    for (i = 0; i < rows; i++) {
        a[i] = (int *)calloc(cols, sizeof(int));
        if (!a[i]) { fprintf(stderr, "alloc2d_int: out of memory\n"); exit(1); }
    }
    return a;
}

static inline void free2d_int(int **a, int rows) {
    int i;
    for (i = 0; i < rows; i++) free(a[i]);
    free(a);
}

static inline long **alloc2d_long(int rows, int cols) {
    int i;
    long **a = (long **)malloc(rows * sizeof(long *));
    if (!a) { fprintf(stderr, "alloc2d_long: out of memory\n"); exit(1); }
    for (i = 0; i < rows; i++) {
        a[i] = (long *)calloc(cols, sizeof(long));
        if (!a[i]) { fprintf(stderr, "alloc2d_long: out of memory\n"); exit(1); }
    }
    return a;
}

static inline void free2d_long(long **a, int rows) {
    int i;
    for (i = 0; i < rows; i++) free(a[i]);
    free(a);
}

/* ─── PairList helpers ───────────────────────────────────────────────────── */
static inline void pl_init(PairList *pl, int cap) {
    pl->p   = (Pair *)malloc(cap * sizeof(Pair));
    if (!pl->p) { fprintf(stderr, "pl_init: out of memory\n"); exit(1); }
    pl->n   = 0;
    pl->cap = cap;
}

static inline void pl_free(PairList *pl) {
    free(pl->p);
    pl->p   = NULL;
    pl->n   = 0;
    pl->cap = 0;
}

static inline void pl_append(PairList *pl, int i, int j) {
    if (pl->n >= pl->cap) {
        pl->cap *= 2;
        pl->p = (Pair *)realloc(pl->p, pl->cap * sizeof(Pair));
        if (!pl->p) { fprintf(stderr, "pl_append: out of memory\n"); exit(1); }
    }
    pl->p[pl->n].i = i;
    pl->p[pl->n].j = j;
    pl->n++;
}

static inline void pl_copy(PairList *dst, const PairList *src) {
    int k;
    if (dst->cap < src->n) {
        dst->p = (Pair *)realloc(dst->p, src->n * sizeof(Pair));
        if (!dst->p) { fprintf(stderr, "pl_copy: out of memory\n"); exit(1); }
        dst->cap = src->n;
    }
    dst->n = src->n;
    for (k = 0; k < src->n; k++) dst->p[k] = src->p[k];
}

/* ─── Function prototypes ────────────────────────────────────────────────── */

/* kepler.c */
void kepler_step(double kc, double dt, State *s0, State *s);
void kepler_step_depth(double kc, double dt, State *s0, State *s, int depth);
int  kepler_step_internal(double kc, double dt, State *s0, State *s);
void kepler_stepxv(double GM, double x[3], double vel[3], double h,
                   double xout[3], double vout[3]);
int  solve_universal_newton(double kc, double r0, double beta,
                            double eta, double zeta, double h,
                            double *X, double *B, double *S2, double *C2);
int  solve_universal_laguerre(double kc, double r0, double beta,
                              double eta, double zeta, double h,
                              double *X, double *B, double *S2, double *C2);
int  solve_universal_newton_bisection(double kc, double r0, double beta,
                                      double eta, double zeta, double h,
                                      double *X, double *B, double *S2, double *C2);
int  solve_universal_parabolic(double kc, double r0, double beta,
                               double eta, double zeta, double h,
                               double *X, double *B, double *S2, double *C2);
int  solve_universal_hyperbolic_newton(double kc, double r0, double minus_beta,
                                       double eta, double zeta, double h, double v2,
                                       double *X, double *B, double *S2, double *C2);
int  solve_universal_hyperbolic_laguerre(double kc, double r0, double minus_beta,
                                         double eta, double zeta, double h, double v2,
                                         double *X, double *B, double *S2, double *C2);
int  solve_universal_hyperbolic_bisection(double kc, double r0, double minus_beta,
                                          double eta, double zeta, double h, double v2,
                                          double *X, double *B, double *S2, double *C2);
double sign_d(double x);
double cubic1(double a, double b, double c);

/* maps.c */
void   map_sun(double **Q, double **v, double h, double *m, int n);
void   interact(double **Q, double **v, double h, double *m, int n, PairList *pairs);
void   adjust_sun(double **v, double *m, int n, double Psum0[3]);
void   calc_m(double **v, double *m, int n, double Psum[3]);
void   calc_levels(double **Q, double **v, double *m, int n,
                   double *rlev, double hglob, int **levc);
void   kep_map(double **Q, double **v, double h,
               double *m, int n, int *indices, int ni);

/* integrator.c */
void   mtr_step(double **Q, double **v, double *m, int n,
                int **levc, double *hlev, double *rlev, double *rlev_sq);
void   driftop(double **Q, double **v, double *m, int n, int lev,
               double *hlev, int **levc, double *rlev, double *rlev_sq,
               PairList *arru);
void   calc_pairs(int lev, int **levc, int n,
                  PairList *arru, PairList *indv, int *kpairs, int *nk);
int    update_levels(double **Q0, double **v0, double **Q, double **v,
                     int **levc, int *levcp,
                     int ind, int hsub, PairList *indv, int lev,
                     PairList *arru, int n, int force_accept);
void   calc_levels_pairs(double **Q, double **v, double *m,
                         double *rlev, double *rlev_sq, double hglob,
                         double *hlev, int **levc,
                         int *levcp, PairList *pairs);

/* coordinates.c */
void   convert_cart(double *m, int n, double **q, double **vel,
                    double **Q, double **P);
void   convert2cart(double *m, int n, double **Q, double **P,
                    double **q, double **p);
void   adjust_cm(double **x, double **v, double *m, int n);

/* physics.c */
void   consqv(double *m, int n, double **Q, double **v,
              double *E, double L[3]);
void   calcorb(double **Q, double **v, double *m,
               int ik, int jk, double *a, double *em);
void   write_snapshot(const char *dir, int isnap, double t,
                      double *m, double **Q, double **v, int n);

/* init.c */
void   initv(double h, double rmax, double rsub, int hsub, int levtot,
             double *rlev, double *hlev, int *steps);
void   in3boddunctwo(double *m, double **q, double **v, int *n);
void   insolardunc(double *m, double **q, double **v, int *n);
void   read_input(const char *fname, SimState *sim, double **q, double **v);

#endif /* MTR_H */
