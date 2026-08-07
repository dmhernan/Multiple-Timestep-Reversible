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

/* --- Snapshot output (Below 2 functions are the ONLY addition to David's code) --- */
static void calcorb_helio(double mu, const double dx[3], const double dv[3],
                          double *a, double *e, double *inc)
{
    double r, v2, rv, hvec[3], hmag, evec[3];
    int k;

    r  = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    v2 = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
    rv = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];

    hvec[0] = dx[1]*dv[2] - dx[2]*dv[1];
    hvec[1] = dx[2]*dv[0] - dx[0]*dv[2];
    hvec[2] = dx[0]*dv[1] - dx[1]*dv[0];
    hmag = sqrt(hvec[0]*hvec[0] + hvec[1]*hvec[1] + hvec[2]*hvec[2]);

    *a = 1.0 / (2.0 / r - v2 / mu);

    for (k = 0; k < 3; k++)
        evec[k] = ((v2 - mu / r) * dx[k] - rv * dv[k]) / mu;
    *e = sqrt(evec[0]*evec[0] + evec[1]*evec[1] + evec[2]*evec[2]);

    *inc = (hmag > 0.0) ? acos(hvec[2] / hmag) : 0.0;
}

/* GPLUM-compatible snapshot (heliocentric, central body NOT included) with (a, e, inc)*/
void write_snapshot(const char *dir, int isnap, double t,
                    double *m, double **Q, double **v, int n)
{
    char fname[512];
    FILE *fp;
    double **x      = alloc2d(3, n);
    double **Pdummy = alloc2d(3, n);
    double **pdummy = alloc2d(3, n);
    double dx[3], dv[3];
    double a, e, inc;
    int i, k;

    convert2cart(m, n, Q, Pdummy, x, pdummy);

    snprintf(fname, sizeof(fname), "%s/snap%06d.dat", dir, isnap);
    fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "write_snapshot: cannot open '%s'\n", fname);
        free2d(x, 3); free2d(Pdummy, 3); free2d(pdummy, 3);
        return;
    }

    fprintf(fp, "%.15e\t%d\t%d", t, n - 1, n - 1);
    for (k = 0; k < 10; k++) fprintf(fp, "\t%.15e", 0.0);
    fprintf(fp, "\n");

    for (i = 1; i < n; i++) {
        for (k = 0; k < 3; k++) {
            dx[k] = x[k][i] - x[k][0];
            dv[k] = v[k][i] - v[k][0];
        }
        calcorb_helio(GNEWT * (m[0] + m[i]), dx, dv, &a, &e, &inc);
        fprintf(fp,
            "%d\t%.15e\t%.15e\t%.15e\t"
            "%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t"
            "%d\t%d\t%.15e\t%.15e\t%.15e\n",
            i - 1, m[i], 0.0, 1.0,
            dx[0], dx[1], dx[2], dv[0], dv[1], dv[2],
            0, 0, a, e, inc);
    }

    fclose(fp);
    free2d(x, 3); free2d(Pdummy, 3); free2d(pdummy, 3);
}
