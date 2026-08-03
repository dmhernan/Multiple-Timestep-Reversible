#include "mtr.h"

void calc_m(double **v, double *m, int n, double Psum[3])
{
    int i, k;
    for (k = 0; k < 3; k++) Psum[k] = 0.0;
    for (i = 0; i < n; i++)
        for (k = 0; k < 3; k++)
            Psum[k] += m[i] * v[k][i];
}

void map_sun(double **Q, double **v, double h, double *m, int n)
{
    double sp[3] = {0.0, 0.0, 0.0};
    int i, k;
    for (i = 1; i < n; i++)
        for (k = 0; k < 3; k++)
            sp[k] += m[i] * v[k][i];
    for (i = 1; i < n; i++)
        for (k = 0; k < 3; k++)
            Q[k][i] += (h / m[0]) * sp[k];
}

void adjust_sun(double **v, double *m, int n, double Psum0[3])
{
    double Psum[3] = {0.0, 0.0, 0.0};
    int i, k;
    for (i = 1; i < n; i++)
        for (k = 0; k < 3; k++)
            Psum[k] += m[i] * v[k][i];
    for (k = 0; k < 3; k++)
        v[k][0] = (Psum0[k] - Psum[k]) / m[0];
}

void interact(double **Q, double **v, double h, double *m, int n, PairList *pairs)
{
    int p, k;
    double qvec[3], qmod2, qmod, fac;
    int ii, ij;
    (void)n;
    if (pairs->n == 0) return;
    for (p = 0; p < pairs->n; p++) {
        ii = pairs->p[p].i;
        ij = pairs->p[p].j;
        qmod2 = 0.0;
        for (k = 0; k < 3; k++) {
            qvec[k] = Q[k][ii] - Q[k][ij];
            qmod2  += qvec[k] * qvec[k];
        }
        qmod = sqrt(qmod2);
        fac  = GNEWT / (qmod2 * qmod);
        for (k = 0; k < 3; k++) {
            v[k][ii] -= h * fac * m[ij] * qvec[k];
            v[k][ij] += h * fac * m[ii] * qvec[k];
        }
    }
}

void kep_map(double **Q, double **v, double h, double *m, int n, int *kpairs, int ni)
{
    double gm = GNEWT * m[0];
    double x[3], vel[3], xout[3], vout[3];
    int i, k, ii;
    (void)n;
    if (ni == 0) return;
    for (i = 0; i < ni; i++) {
        ii = kpairs[i];
        for (k = 0; k < 3; k++) {
            x[k]   = Q[k][ii];
            vel[k] = v[k][ii];
        }
        kepler_stepxv(gm, x, vel, h, xout, vout);
        for (k = 0; k < 3; k++) {
            Q[k][ii] = xout[k];
            v[k][ii] = vout[k];
        }
    }
}

void calc_levels(double **Q, double **v, double *m, int n,
                 double *rlev, double hglob, int **levc)
{
    int j, k, lev, i;
    double Qjk[3], qjk_sq, ratio2;
    int levtot = LEVMAX;
    double hglob2 = hglob * hglob;
    double rhill2 = G_RHILL * G_RHILL;
    (void)v;

    for (j = 0; j < n; j++)
        for (k = 0; k < n; k++)
            levc[j][k] = 0;

    for (j = 1; j < n; j++) {
        for (k = j+1; k < n; k++) {
            qjk_sq = 0.0;
            for (i = 0; i < 3; i++) {
                Qjk[i] = Q[i][j] - Q[i][k];
                qjk_sq += Qjk[i] * Qjk[i];
            }

            if (G_FLAG == 0) {
                double msum = GNEWT * (m[j] + m[k]) * hglob2;
                ratio2 = qjk_sq * sqrt(qjk_sq) / msum;
            } else {
                /* flag==1 and flag==2: use Q_ij/r_Hill for initial assignment */
                ratio2 = qjk_sq / rhill2;
            }

            if (ratio2 > rlev[0]*rlev[0]) {
                lev = 1;
            } else {
                lev = 0;
                for (i = 0; i < levtot - 1; i++) {
                    if (ratio2 <= rlev[i]*rlev[i] && ratio2 > rlev[i+1]*rlev[i+1]) {
                        lev = i + 1;
                        break;
                    }
                }
            }

            if (lev > 0) {
                levc[j][k] = lev;
                levc[k][j] = lev;
            }
        }
    }
}

/* Targeted version: only recompute levels for pairs in the given PairList.
   rlev_sq[i] = rlev[i]*rlev[i] precomputed by caller.
   For G_FLAG==0: ratio^2 = qjk^3/(G*M*hglob^2), needs one sqrt (qjk).
   For G_FLAG==1: ratio^2 = qjk_sq/G_RHILL^2, needs zero sqrts.
   For G_FLAG==2: ratio^2 = qmin^2/G_RHILL^2, using hlev[levc[j][k]-1] per pair. */
void calc_levels_pairs(double **Q, double **v, double *m,
                       double *rlev, double *rlev_sq, double hglob,
                       double *hlev, int **levc,
                       int *levcp, PairList *pairs)
{
    int p, i, lev_idx;
    double Qjk[3], Vjk[3], qjk_sq, ratio2;
    int levtot = LEVMAX;
    double hglob2 = hglob * hglob;
    double rhill2 = G_RHILL * G_RHILL;
    (void)rlev;

    for (p = 0; p < pairs->n; p++) {
        int j = pairs->p[p].i;
        int k = pairs->p[p].j;

        qjk_sq = 0.0;
        for (i = 0; i < 3; i++) {
            Qjk[i] = Q[i][j] - Q[i][k];
            qjk_sq += Qjk[i] * Qjk[i];
        }

        if (G_FLAG == 0) {
            double msum = GNEWT * (m[j] + m[k]) * hglob2;
            ratio2 = qjk_sq * sqrt(qjk_sq) / msum;
        } else if (G_FLAG == 1) {
            ratio2 = qjk_sq / rhill2;
        } else {
            /* flag==2: Q_min/r_Hill, using this pair's current level timestep */
            double hi, QdotV, V2, qmin2, d, tmin;
            int cur_lev = (levc != NULL && levc[j][k] >= 1) ? levc[j][k] : 1;
            hi = (hlev != NULL) ? hlev[cur_lev - 1] : hglob;

            for (i = 0; i < 3; i++)
                Vjk[i] = v[i][j] - v[i][k];

            QdotV = 0.0; V2 = 0.0;
            for (i = 0; i < 3; i++) {
                QdotV += Qjk[i] * Vjk[i];
                V2    += Vjk[i] * Vjk[i];
            }

            if (V2 == 0.0 || QdotV == 0.0) {
                qmin2 = qjk_sq;
            } else {
                d    = (QdotV > 0.0) ? -1.0 : 1.0;
                tmin = -d * QdotV / V2;
                if (tmin < hi / 2.0) {
                    qmin2 = qjk_sq - QdotV * QdotV / V2;
                } else {
                    qmin2 = qjk_sq + hi * d * QdotV + hi * hi * V2 / 4.0;
                }
                if (qmin2 < 0.0) qmin2 = 0.0;
            }
            ratio2 = qmin2 / rhill2;
        }

        if (ratio2 > rlev_sq[0]) {
            lev_idx = 1;
        } else {
            lev_idx = 0;
            for (i = 0; i < levtot - 1; i++) {
                if (ratio2 <= rlev_sq[i] && ratio2 > rlev_sq[i+1]) {
                    lev_idx = i + 1;
                    break;
                }
            }
        }

        levcp[p] = lev_idx;
    }
}
