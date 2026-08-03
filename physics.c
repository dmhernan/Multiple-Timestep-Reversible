#include "mtr.h"

static void cross3(double a[3], double b[3], double c[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

void consqv(double *m, int n, double **Q, double **v, double *E, double L[3])
{
    double **x       = alloc2d(3, n);
    double **Pdummy  = alloc2d(3, n);
    double **pdummy  = alloc2d(3, n);
    double T, V;
    double xij[3], rij;
    int i, j, k;

    convert2cart(m, n, Q, Pdummy, x, pdummy);

    T = 0.0; V = 0.0;
    for (k = 0; k < 3; k++) L[k] = 0.0;

    for (i = 0; i < n; i++) {
        double vi2 = 0.0;
        double xi[3], vi[3], lci[3];
        for (k = 0; k < 3; k++) {
            xi[k] = x[k][i];
            vi[k] = v[k][i];
            vi2  += vi[k] * vi[k];
        }
        T += 0.5 * m[i] * vi2;
        cross3(xi, vi, lci);
        for (k = 0; k < 3; k++) L[k] += m[i] * lci[k];
        for (j = i+1; j < n; j++) {
            rij = 0.0;
            for (k = 0; k < 3; k++) {
                xij[k] = x[k][i] - x[k][j];
                rij   += xij[k] * xij[k];
            }
            rij = sqrt(rij);
            V  -= GNEWT * m[i] * m[j] / rij;
        }
    }
    *E = T + V;

    free2d(x, 3);
    free2d(Pdummy, 3);
    free2d(pdummy, 3);
}

void calcorb(double **Q, double **v, double *m, int ik, int jk,
             double *a, double *em)
{
    double qr[3], vr[3];
    double qrm, vrm2, M, mudan, mu;
    double h[3], pr[3], e[3], hc[3], tmp[3];
    int k;

    M     = m[ik] + m[jk];
    mu    = m[ik] * m[jk] / M;
    mudan = GNEWT * M;

    qrm = 0.0; vrm2 = 0.0;
    for (k = 0; k < 3; k++) {
        qr[k]  = Q[k][ik] - Q[k][jk];
        vr[k]  = v[k][ik] - v[k][jk];
        qrm   += qr[k] * qr[k];
        vrm2  += vr[k] * vr[k];
    }
    qrm = sqrt(qrm);

    for (k = 0; k < 3; k++) pr[k] = mu * vr[k];
    cross3(qr, vr, hc);
    for (k = 0; k < 3; k++) h[k] = mu * hc[k];
    cross3(pr, h, tmp);
    for (k = 0; k < 3; k++)
        e[k] = (tmp[k] - GNEWT * mu * mu * M * qr[k] / qrm)
               / (GNEWT * mu * mu * M);

    *em = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
    *a  = 1.0 / (2.0 / qrm - vrm2 / mudan);
}
