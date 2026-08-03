#include "mtr.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ═══════════════════════════════════════════════════════════════════════════
   kepler.c
   Kepler drift operator — universal variable formulation.
   Root-finding solvers follow Wisdom & Hernandez (2015).
   The f/g application uses the fhat/gdothat form from the MATLAB source
   to avoid catastrophic cancellation when the timestep is small.
   ═══════════════════════════════════════════════════════════════════════════ */

/* ─── Helpers ────────────────────────────────────────────────────────────── */

double sign_d(double x)
{
    if (x > 0.0) return  1.0;
    if (x < 0.0) return -1.0;
    return 0.0;
}

/*  Real root of  t^3 + (a/3)t^2 + (b/3)t + c = 0  via Cardano.
    Exits if three real roots exist (not expected in orbital mechanics). */
double cubic1(double a, double b, double c)
{
    double Q = (a*a - 3.0*b) / 9.0;
    double R = (2.0*a*a*a - 9.0*a*b + 27.0*c) / 54.0;
    double A, B;

    if (R*R < Q*Q*Q) {
        fprintf(stderr, "cubic1: three real roots, exiting.\n");
        exit(-1);
    }
    A = -sign_d(R) * pow(fabs(R) + sqrt(R*R - Q*Q*Q), 1.0/3.0);
    B = (A == 0.0) ? 0.0 : Q / A;
    return (A + B) - a / 3.0;
}

/* ─── Elliptic solvers ───────────────────────────────────────────────────── */

int solve_universal_newton(double kc, double r0, double beta,
                           double eta, double zeta, double h,
                           double *X, double *B, double *S2, double *C2)
{
    (void)kc;
    double b    = sqrt(beta);
    double xnew = *X;
    double err  = 1.e-9 * fabs(xnew);
    double x, arg, s2, c2, g1, g2, g3, g;
    int count = 0;

    do {
        x    = xnew;
        arg  = b * x / 2.0;
        s2   = sin(arg);
        c2   = cos(arg);
        g1   = 2.0 * s2 * c2 / b;
        g2   = 2.0 * s2 * s2 / beta;
        g3   = (x - g1) / beta;
        g    = eta * g1 + zeta * g2;
        xnew = ((x * g - eta * g2 - zeta * g3) + h) / (r0 + g);
        if (count++ > 10) return FAILURE;
    } while (fabs(x - xnew) > err);

    x    = xnew;
    arg  = b * x / 2.0;
    *X   = x;
    *B   = b;
    *S2  = sin(arg);
    *C2  = cos(arg);
    return SUCCESS;
}

int solve_universal_laguerre(double kc, double r0, double beta,
                             double eta, double zeta, double h,
                             double *X, double *B, double *S2, double *C2)
{
    (void)kc;
    double b    = sqrt(beta);
    double xnew = *X;
    double x, arg, s2, c2, g1, g2, g3, g0, f, fp, fpp, dx;
    int count = 0;
    const double c5 = 5.0, c16 = 16.0, c20 = 20.0;

    do {
        x    = xnew;
        arg  = b * x / 2.0;
        s2   = sin(arg);
        c2   = cos(arg);
        g1   = 2.0 * s2 * c2 / b;
        g2   = 2.0 * s2 * s2 / beta;
        g3   = (x - g1) / beta;
        f    = r0 * x + eta * g2 + zeta * g3 - h;
        fp   = r0 + eta * g1 + zeta * g2;
        g0   = 1.0 - beta * g2;
        fpp  = eta * g0 + zeta * g1;
        dx   = -c5 * f / (fp + sqrt(fabs(c16 * fp * fp - c20 * f * fpp)));
        xnew = x + dx;
        if (count++ > 20) return FAILURE;
    } while (fabs(dx) > 1.e-9 * fabs(xnew));

    x    = xnew;
    arg  = b * x / 2.0;
    *X   = x;
    *B   = b;
    *S2  = sin(arg);
    *C2  = cos(arg);
    return SUCCESS;
}

int solve_universal_newton_bisection(double kc, double r0, double beta,
                                     double eta, double zeta, double h,
                                     double *X, double *B, double *S2, double *C2)
{
    (void)kc;
    double b    = sqrt(beta);
    double xnew = *X;
    double err  = 1.e-9 * fabs(xnew);
    double x, arg, s2, c2, g1, g2, g3, g, f;
    double X_min, X_max, invperiod, X_per_period;
    int count = 0;

    /* Bisection bounds: Rein et al. */
    invperiod    = b * beta / (2.0 * M_PI * kc);
    X_per_period = 2.0 * M_PI / b;
    X_min = X_per_period * floor(h * invperiod);
    X_max = X_min + X_per_period;
    xnew  = (X_max + X_min) / 2.0;

    do {
        x   = xnew;
        arg = b * x / 2.0;
        s2  = sin(arg);
        c2  = cos(arg);
        g1  = 2.0 * s2 * c2 / b;
        g2  = 2.0 * s2 * s2 / beta;
        g3  = (x - g1) / beta;
        f   = r0 * x + eta * g2 + zeta * g3 - h;
        if (f >= 0.0) X_max = x;
        else          X_min = x;
        xnew = (X_max + X_min) / 2.0;
        if (count++ > 100) return FAILURE;
    } while (fabs(x - xnew) > err);

    /* Final Newton polish */
    x   = xnew;
    arg = b * x / 2.0;
    s2  = sin(arg);
    c2  = cos(arg);
    g1  = 2.0 * s2 * c2 / b;
    g2  = 2.0 * s2 * s2 / beta;
    g3  = (x - g1) / beta;
    g   = eta * g1 + zeta * g2;
    x   = ((x * g - eta * g2 - zeta * g3) + h) / (r0 + g);

    arg  = b * x / 2.0;
    *X   = x;
    *B   = b;
    *S2  = sin(arg);
    *C2  = cos(arg);
    return SUCCESS;
}

/* ─── Parabolic solver ───────────────────────────────────────────────────── */

int solve_universal_parabolic(double kc, double r0, double beta,
                              double eta, double zeta, double h,
                              double *X, double *B, double *S2, double *C2)
{
    (void)kc;
    (void)beta;
    double x = cubic1(3.0 * eta / zeta, 6.0 * r0 / zeta, -6.0 * h / zeta);
    *X  = x;
    *B  = 0.0;
    *S2 = 0.0;
    *C2 = 1.0;
    return SUCCESS;
}

/* ─── Hyperbolic solvers ─────────────────────────────────────────────────── */

int solve_universal_hyperbolic_newton(double kc, double r0, double minus_beta,
                                      double eta, double zeta, double h, double v2,
                                      double *X, double *B, double *S2, double *C2)
{
    (void)v2;
    (void)kc;
    double b    = sqrt(minus_beta);
    double xnew = *X;
    double x, arg, s2, c2, g1, g2, g3, g;
    int count = 0;

    do {
        x   = xnew;
        arg = b * x / 2.0;
        if (fabs(arg) > 50.0) return FAILURE;
        s2  = sinh(arg);
        c2  = cosh(arg);
        g1  = 2.0 * s2 * c2 / b;
        g2  = 2.0 * s2 * s2 / minus_beta;
        g3  = -(x - g1) / minus_beta;
        g   = eta * g1 + zeta * g2;
        xnew = (x * g - eta * g2 - zeta * g3 + h) / (r0 + g);
        if (count++ > 10) return FAILURE;
    } while (fabs(x - xnew) > 1.e-9 * fabs(xnew));

    /* Final Newton polish */
    x   = xnew;
    arg = b * x / 2.0;
    s2  = sinh(arg);
    c2  = cosh(arg);
    g1  = 2.0 * s2 * c2 / b;
    g2  = 2.0 * s2 * s2 / minus_beta;
    g3  = -(x - g1) / minus_beta;
    g   = eta * g1 + zeta * g2;
    x   = (x * g - eta * g2 - zeta * g3 + h) / (r0 + g);

    arg  = b * x / 2.0;
    *X   = x;
    *B   = b;
    *S2  = sinh(arg);
    *C2  = cosh(arg);
    return SUCCESS;
}

int solve_universal_hyperbolic_laguerre(double kc, double r0, double minus_beta,
                                        double eta, double zeta, double h, double v2,
                                        double *X, double *B, double *S2, double *C2)
{
    (void)v2;
    (void)kc;
    double b    = sqrt(minus_beta);
    double xnew = *X;
    double x, arg, s2, c2, g0, g1, g2, g3, f, fp, fpp, dx, den;
    int count = 0;
    const double c5 = 5.0, c16 = 16.0, c20 = 20.0;

    do {
        x   = xnew;
        arg = b * x / 2.0;
        if (fabs(arg) > 50.0) return FAILURE;
        s2  = sinh(arg);
        c2  = cosh(arg);
        g1  = 2.0 * s2 * c2 / b;
        g2  = 2.0 * s2 * s2 / minus_beta;
        g3  = -(x - g1) / minus_beta;
        f   = r0 * x + eta * g2 + zeta * g3 - h;
        fp  = r0 + eta * g1 + zeta * g2;
        g0  = 1.0 + minus_beta * g2;
        fpp = eta * g0 + zeta * g1;
        den = fp + sqrt(fabs(c16 * fp * fp - c20 * f * fpp));
        if (den == 0.0) return FAILURE;
        dx   = -c5 * f / den;
        xnew = x + dx;
        if (count++ > 20) return FAILURE;
    } while (fabs(x - xnew) > 1.e-9 * fabs(xnew));

    /* Final Newton polish */
    x   = xnew;
    arg = b * x / 2.0;
    if (fabs(arg) > 200.0) return FAILURE;
    s2  = sinh(arg);
    c2  = cosh(arg);
    g1  = 2.0 * s2 * c2 / b;
    g2  = 2.0 * s2 * s2 / minus_beta;
    g3  = -(x - g1) / minus_beta;
    {
        double g = eta * g1 + zeta * g2;
        x = (x * g - eta * g2 - zeta * g3 + h) / (r0 + g);
    }

    arg  = b * x / 2.0;
    *X   = x;
    *B   = b;
    *S2  = sinh(arg);
    *C2  = cosh(arg);
    return SUCCESS;
}

int solve_universal_hyperbolic_bisection(double kc, double r0, double minus_beta,
                                         double eta, double zeta, double h, double v2,
                                         double *X, double *B, double *S2, double *C2)
{
    (void)v2;
    (void)kc;
    double b    = sqrt(minus_beta);
    double xnew = *X;
    double x, arg, s2, c2, g1, g2, g3, g, f, fmin, fmax;
    double X_min, X_max;
    double h2, q, vq;
    int count = 0;

    /* Bounds from Wisdom & Hernandez */
    h2    = r0 * r0 * v2 - eta * eta;
    q     = h2 / kc / (1.0 + sqrt(1.0 + h2 * minus_beta / (kc * kc)));
    vq    = sqrt(h2) / q;
    X_min = 0.99 / (vq + r0 / h);   /* slight kludge per W&H */
    X_max = 10.0 * xnew;

    if ((X_min - xnew) * (X_max - xnew) > 0.0) return FAILURE;

    /* Evaluate at bounds */
    x   = X_min;
    arg = b * x / 2.0;
    if (fabs(arg) > 50.0) return FAILURE;
    s2   = sinh(arg); c2 = cosh(arg);
    g1   = 2.0 * s2 * c2 / b;
    g2   = 2.0 * s2 * s2 / minus_beta;
    g3   = -(x - g1) / minus_beta;
    fmin = r0 * x + eta * g2 + zeta * g3 - h;

    x   = X_max;
    arg = b * x / 2.0;
    if (fabs(arg) > 50.0) { x = 50.0 / (b / 2.0); arg = 50.0; }
    s2   = sinh(arg); c2 = cosh(arg);
    g1   = 2.0 * s2 * c2 / b;
    g2   = 2.0 * s2 * s2 / minus_beta;
    g3   = -(x - g1) / minus_beta;
    fmax = r0 * x + eta * g2 + zeta * g3 - h;

    if (fmin * fmax > 0.0) return FAILURE;

    do {
        x   = xnew;
        arg = b * x / 2.0;
        if (fabs(arg) > 50.0) return FAILURE;
        s2  = sinh(arg); c2 = cosh(arg);
        g1  = 2.0 * s2 * c2 / b;
        g2  = 2.0 * s2 * s2 / minus_beta;
        g3  = -(x - g1) / minus_beta;
        f   = r0 * x + eta * g2 + zeta * g3 - h;
        if (f >= 0.0) X_max = x;
        else          X_min = x;
        xnew = (X_max + X_min) / 2.0;
        if (count++ > 100) return FAILURE;
    } while (fabs(x - xnew) > 1.e-12 * fabs(xnew));

    /* Final Newton polish */
    x   = xnew;
    arg = b * x / 2.0;
    if (fabs(arg) > 50.0) return FAILURE;
    s2  = sinh(arg); c2 = cosh(arg);
    g1  = 2.0 * s2 * c2 / b;
    g2  = 2.0 * s2 * s2 / minus_beta;
    g3  = -(x - g1) / minus_beta;
    g   = eta * g1 + zeta * g2;
    x   = (x * g - eta * g2 - zeta * g3 + h) / (r0 + g);

    arg  = b * x / 2.0;
    *X   = x;
    *B   = b;
    *S2  = sinh(arg);
    *C2  = cosh(arg);
    return SUCCESS;
}

/* ─── Initial guess ──────────────────────────────────────────────────────── */

static double new_guess(double r0, double eta, double zeta, double dt)
{
    if (zeta != 0.0) {
        return cubic1(3.0 * eta / zeta, 6.0 * r0 / zeta, -6.0 * dt / zeta);
    } else if (eta != 0.0) {
        double reta = r0 / eta;
        double disc = reta * reta + 8.0 * dt / eta;
        return (disc >= 0.0) ? (-reta + sqrt(disc)) : dt / r0;
    } else {
        return dt / r0;
    }
}

/* ─── Internal Kepler step ───────────────────────────────────────────────── */

int kepler_step_internal(double kc, double dt, State *s0, State *s)
{
    double r0, v2, eta, beta, zeta, b;
    double s2, c2, G1, G2, r, a, c, ca, bsa;
    double fhat, gdothat, g;
    double x, x0;
    int flag;

    r0   = sqrt(s0->x*s0->x + s0->y*s0->y + s0->z*s0->z);
    v2   = s0->xd*s0->xd + s0->yd*s0->yd + s0->zd*s0->zd;
    eta  = s0->x*s0->xd + s0->y*s0->yd + s0->z*s0->zd;
    beta = 2.0 * kc / r0 - v2;
    zeta = kc - beta * r0;

    if (beta < 0.0) {
        /* Hyperbolic */
        x0   = new_guess(r0, eta, zeta, dt);
        x    = x0;
        flag = solve_universal_hyperbolic_newton(kc, r0, -beta, eta, zeta, dt, v2,
                                                 &x, &b, &s2, &c2);
        if (flag == FAILURE) {
            x    = x0;
            flag = solve_universal_hyperbolic_laguerre(kc, r0, -beta, eta, zeta, dt, v2,
                                                       &x, &b, &s2, &c2);
        }
        if (flag == FAILURE) {
            x    = x0;
            flag = solve_universal_hyperbolic_bisection(kc, r0, -beta, eta, zeta, dt, v2,
                                                        &x, &b, &s2, &c2);
        }
        if (flag == FAILURE) return FAILURE;

        a   = kc / (-beta);
        G1  = 2.0 * s2 * c2 / b;
        c   = 2.0 * s2 * s2;
        G2  = c / (-beta);
        ca  = c * a;
        r   = r0 + eta * G1 + zeta * G2;
        bsa = (a / r) * (b / r0) * 2.0 * s2 * c2;

    } else if (beta > 0.0) {
        /* Elliptic */
        x0 = dt / r0;
        {
            double ff = zeta * x0*x0*x0 + 3.0 * eta * x0*x0;
            double fp = 3.0 * zeta * x0*x0 + 6.0 * eta * x0 + 6.0 * r0;
            x0 = x0 - ff / fp;
        }
        x    = x0;
        flag = solve_universal_newton(kc, r0, beta, eta, zeta, dt,
                                      &x, &b, &s2, &c2);
        if (flag == FAILURE) {
            x    = x0;
            flag = solve_universal_laguerre(kc, r0, beta, eta, zeta, dt,
                                            &x, &b, &s2, &c2);
        }
        if (flag == FAILURE) {
            x    = x0;
            flag = solve_universal_newton_bisection(kc, r0, beta, eta, zeta, dt,
                                                    &x, &b, &s2, &c2);
        }
        if (flag == FAILURE) return FAILURE;

        a   = kc / beta;
        G1  = 2.0 * s2 * c2 / b;
        c   = 2.0 * s2 * s2;
        G2  = c / beta;
        ca  = c * a;
        r   = r0 + eta * G1 + zeta * G2;
        bsa = (a / r) * (b / r0) * 2.0 * s2 * c2;

    } else {
        /* Parabolic */
        x    = dt / r0;
        flag = solve_universal_parabolic(kc, r0, beta, eta, zeta, dt,
                                         &x, &b, &s2, &c2);
        if (flag == FAILURE) { fprintf(stderr, "parabolic solver failed\n"); exit(-1); }

        G1  = x;
        G2  = x * x / 2.0;
        ca  = kc * G2;
        r   = r0 + eta * G1 + zeta * G2;
        bsa = kc * x / (r * r0);
    }

    /* ── Apply f/g map using fhat = f-1, gdothat = gdot-1
          to avoid cancellation when timestep is small (MATLAB version) ── */
    fhat    = -(ca / r0);          /*  = f - 1          */
    g       = eta * G2 + r0 * G1;
    gdothat = -(ca / r);           /*  = gdot - 1       */

    s->x  = s0->x  + fhat * s0->x  + g * s0->xd;
    s->y  = s0->y  + fhat * s0->y  + g * s0->yd;
    s->z  = s0->z  + fhat * s0->z  + g * s0->zd;
    s->xd = s0->xd - bsa  * s0->x  + gdothat * s0->xd;
    s->yd = s0->yd - bsa  * s0->y  + gdothat * s0->yd;
    s->zd = s0->zd - bsa  * s0->z  + gdothat * s0->zd;

    return SUCCESS;
}

/* ─── Recursive depth-halving wrapper ───────────────────────────────────── */

void kepler_step_depth(double kc, double dt, State *s0, State *s, int depth)
{
    State ss;

    if (depth > 30) {
        fprintf(stderr, "kepler_step_depth: depth exceeded\n");
        exit(-1);
    }

    if (kepler_step_internal(kc, dt, s0, s) == SUCCESS) return;

    /* Subdivide into 4 substeps of dt/4 */
    kepler_step_depth(kc, dt / 4.0, s0,  &ss, depth + 1);
    kepler_step_depth(kc, dt / 4.0, &ss, s,   depth + 1);
    kepler_step_depth(kc, dt / 4.0, s,   &ss, depth + 1);
    kepler_step_depth(kc, dt / 4.0, &ss, s,   depth + 1);
}

void kepler_step(double kc, double dt, State *s0, State *s)
{
    kepler_step_depth(kc, dt, s0, s, 0);
}

/* ─── Array interface (used by kep_map in maps.c) ─────────────────────── */

/*  Takes position x[3] and velocity vel[3] (column vectors as in MATLAB),
    advances by timestep h under central force GM, returns xout[3], vout[3]. */
void kepler_stepxv(double GM, double x[3], double vel[3], double h,
                   double xout[3], double vout[3])
{
    State s0, s;

    s0.x  = x[0];   s0.y  = x[1];   s0.z  = x[2];
    s0.xd = vel[0]; s0.yd = vel[1]; s0.zd = vel[2];

    kepler_step(GM, h, &s0, &s);

    xout[0] = s.x;  xout[1] = s.y;  xout[2] = s.z;
    vout[0] = s.xd; vout[1] = s.yd; vout[2] = s.zd;
}
