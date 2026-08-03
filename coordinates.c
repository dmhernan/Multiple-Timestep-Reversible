#include "mtr.h"

void convert_cart(double *m, int n, double **q, double **vel,
                  double **Q, double **P)
{
    double M = 0.0;
    double *a = (double *)malloc(n * sizeof(double));
    double **p = alloc2d(3, n);
    int i, k;

    for (i = 0; i < n; i++) M += m[i];
    for (i = 0; i < n; i++) a[i] = m[i] / M;

    for (i = 0; i < n; i++)
        for (k = 0; k < 3; k++)
            p[k][i] = m[i] * vel[k][i];

    for (k = 0; k < 3; k++) {
        Q[k][0] = 0.0; P[k][0] = 0.0;
        for (i = 0; i < n; i++) {
            Q[k][0] += a[i] * q[k][i];
            P[k][0] += p[k][i];
        }
    }

    for (i = 1; i < n; i++)
        for (k = 0; k < 3; k++) {
            P[k][i] = p[k][i] - a[i] * P[k][0];
            Q[k][i] = q[k][i] - q[k][0];
        }

    free(a);
    free2d(p, 3);
}

void convert2cart(double *m, int n, double **Q, double **P,
                  double **q, double **p)
{
    double M = 0.0;
    double *a = (double *)malloc(n * sizeof(double));
    double s1[3], s2[3];
    int i, k;

    for (i = 0; i < n; i++) M += m[i];
    for (i = 0; i < n; i++) a[i] = m[i] / M;

    for (k = 0; k < 3; k++) {
        s1[k] = 0.0; s2[k] = 0.0;
        for (i = 1; i < n; i++) {
            s1[k] += a[i] * Q[k][i];
            s2[k] += P[k][i];
        }
    }

    for (k = 0; k < 3; k++) {
        p[k][0] = a[0] * P[k][0] - s2[k];
        q[k][0] = Q[k][0] - s1[k];
    }

    for (i = 1; i < n; i++)
        for (k = 0; k < 3; k++) {
            p[k][i] = a[i] * P[k][0] + P[k][i];
            q[k][i] = Q[k][0] + Q[k][i] - s1[k];
        }

    free(a);
}

void adjust_cm(double **x, double **v, double *m, int n)
{
    double ms = 0.0;
    double vcm[3] = {0.0, 0.0, 0.0};
    double xcm[3] = {0.0, 0.0, 0.0};
    int i, k;

    for (i = 0; i < n; i++) ms += m[i];
    for (i = 0; i < n; i++)
        for (k = 0; k < 3; k++) {
            vcm[k] += m[i] * v[k][i];
            xcm[k] += m[i] * x[k][i];
        }
    for (k = 0; k < 3; k++) { vcm[k] /= ms; xcm[k] /= ms; }
    for (i = 0; i < n; i++)
        for (k = 0; k < 3; k++) {
            x[k][i] -= xcm[k];
            v[k][i] -= vcm[k];
        }
}
