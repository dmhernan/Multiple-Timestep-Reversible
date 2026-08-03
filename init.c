#include "mtr.h"

void initv(double h, double rmax, double rsub, int hsub, int levtot,
           double *rlev, double *hlev, int *steps)
{
    double rmin = 0.0;
    double delr = rmax - rmin;
    int i;

    rlev[0]  = rmax;
    hlev[0]  = h;
    steps[0] = 1;

    for (i = 1; i < levtot; i++) {
        rlev[i]  = delr / pow(rsub, (double)i);
        hlev[i]  = h / pow((double)hsub, (double)i);
        steps[i] = steps[i-1] * hsub;
    }
}

void in3boddunctwo(double *m, double **q, double **v, int *n)
{
    double a0  = 0.0125, a1  = 1.0,  e0  = 0.6, e1  = 0.0;
    double a2  = 0.0130, a1p = 3.0,  e2  = 0.2, e1p = 0.0;
    double m0  = 1e-3;
    double mt0, mt1;
    double x0, x1, x1p, x2, v0, v1, v1p, v2;
    double fa0, fb0, fa1, fb1;
    int k, i;
    double **qtmp = alloc2d(3, 5);
    double **vtmp = alloc2d(3, 5);

    *n = 5;
    m[0] = 1.0; m[1] = m0; m[2] = m0; m[3] = m0; m[4] = m0;

    mt0 = m[1] + m[2];
    mt1 = m[0] + mt0;

    x0  = a0  * (1.0 + e0);
    x1  = a1  * (1.0 + e1);
    x1p = a1p * (1.0 + e1p);
    x2  = a2  * (1.0 + e2);

    v0  = sqrt(GNEWT * mt0 / a0  * (1.0 - e0)  / (1.0 + e0));
    v1  = sqrt(GNEWT * mt1 / a1  * (1.0 - e1)  / (1.0 + e1));
    v1p = sqrt(GNEWT * mt1 / a1p * (1.0 - e1p) / (1.0 + e1p));
    v2  = sqrt(GNEWT * mt0 / a2  * (1.0 - e2)  / (1.0 + e2));

    fa0 = m[1] / mt0; fb0 = m[2] / mt0;
    fa1 = mt0  / mt1; fb1 = m[0] / mt1;

    for (k = 0; k < 3; k++)
        for (i = 0; i < 5; i++) { qtmp[k][i] = 0.0; vtmp[k][i] = 0.0; }

    qtmp[0][0] =  fa1 * x1;
    qtmp[0][1] = -fb1 * x1 - fb0 * x0;
    qtmp[0][2] = -fb1 * x1 + fa0 * x0;
    qtmp[0][3] = -fb1 * x1p - fb0 * x2;
    qtmp[0][4] = -fb1 * x1p + fa0 * x2;

    vtmp[1][0] =  fa1 * v1;
    vtmp[1][1] = -fb1 * v1 - fb0 * v0;
    vtmp[1][2] = -fb1 * v1 + fa0 * v0;
    vtmp[1][3] = -fb1 * v1p - fb0 * v2;
    vtmp[1][4] = -fb1 * v1p + fa0 * v2;

    adjust_cm(qtmp, vtmp, m, *n);

    for (k = 0; k < 3; k++)
        for (i = 0; i < *n; i++) {
            q[k][i] = qtmp[k][i];
            v[k][i] = vtmp[k][i];
        }

    free2d(qtmp, 3);
    free2d(vtmp, 3);
}

void insolardunc(double *m, double **q, double **v, int *n)
{
    double fac = 50.0;
    int i, k;
    double vcm[3], xcm[3], ms;

    *n = 5;

    m[0] = 1.00000597682;
    m[1] = fac * 0.000954786104043;
    m[2] = fac * 0.000285583733151;
    m[3] = fac * 0.0000437273164546;
    m[4] = fac * 0.0000517759138449;

    q[0][0]=0.0;         q[1][0]=0.0;          q[2][0]=0.0;
    q[0][1]=-3.5023653;  q[1][1]=-3.8169847;   q[2][1]=-1.5507963;
    q[0][2]= 9.0755314;  q[1][2]=-3.0458353;   q[2][2]=-1.6483708;
    q[0][3]= 8.3101420;  q[1][3]=-16.2901086;  q[2][3]=-7.2521278;
    q[0][4]=11.4707666;  q[1][4]=-25.7294829;  q[2][4]=-10.8169456;

    v[0][0]=0.0;         v[1][0]=0.0;          v[2][0]=0.0;
    v[0][1]= 0.00565429; v[1][1]=-0.00412490;  v[2][1]=-0.00190589;
    v[0][2]= 0.00168318; v[1][2]= 0.00483525;  v[2][2]= 0.00192462;
    v[0][3]= 0.00354178; v[1][3]= 0.00137102;  v[2][3]= 0.00055029;
    v[0][4]= 0.00288930; v[1][4]= 0.00114527;  v[2][4]= 0.00039677;

    ms = 0.0;
    for (k = 0; k < 3; k++) { vcm[k] = 0.0; xcm[k] = 0.0; }
    for (i = 0; i < *n; i++) {
        ms += m[i];
        for (k = 0; k < 3; k++) {
            vcm[k] += m[i] * v[k][i];
            xcm[k] += m[i] * q[k][i];
        }
    }
    for (k = 0; k < 3; k++) { vcm[k] /= ms; xcm[k] /= ms; }
    for (i = 0; i < *n; i++)
        for (k = 0; k < 3; k++) {
            v[k][i] -= vcm[k];
            q[k][i] -= xcm[k];
        }
}

void read_input(const char *fname, SimState *sim, double **q, double **v)
{
    FILE *fp;
    char line[512], key[64];
    double val;
    int i;
    double rmax = 30.0, rsub = 2.0;
    int levtot = 30;
    int steps[LEVMAX];

    fp = fopen(fname, "r");
    if (!fp) { fprintf(stderr, "read_input: cannot open '%s'\n", fname); exit(-1); }

    sim->n = 0; sim->h = 1.0; sim->tmax = 365250.0;
    sim->output_interval = 1; sim->flag = 0;
    sim->hsub = 3; sim->rhill = 0.0;

    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        if (sscanf(line, "%s %lf", key, &val) == 2) {
            if      (strcmp(key, "n_bodies")        == 0) sim->n               = (int)val;
            else if (strcmp(key, "timestep")        == 0) sim->h               = val;
            else if (strcmp(key, "tmax")            == 0) sim->tmax            = val;
            else if (strcmp(key, "output_interval") == 0) sim->output_interval = (int)val;
            else if (strcmp(key, "rmax")            == 0) rmax                 = val;
            else if (strcmp(key, "rsub")            == 0) rsub                 = val;
            else if (strcmp(key, "hsub")            == 0) sim->hsub            = (int)val;
            else if (strcmp(key, "levtot")          == 0) levtot               = (int)val;
            else if (strcmp(key, "flag")            == 0) sim->flag            = (int)val;
            else if (strcmp(key, "rhill")           == 0) sim->rhill           = val;
        }
    }

    rewind(fp);
    i = 0;
    while (fgets(line, sizeof(line), fp) && i < sim->n) {
        if (line[0] == '#' || line[0] == '\n') continue;
        if (strchr("0123456789-+.", line[0])) {
            if (sscanf(line, "%lf %lf %lf %lf %lf %lf %lf",
                       &sim->m[i],
                       &q[0][i], &q[1][i], &q[2][i],
                       &v[0][i], &v[1][i], &v[2][i]) == 7)
                i++;
        }
    }
    fclose(fp);

    if (i != sim->n) {
        fprintf(stderr, "read_input: expected %d bodies, got %d\n", sim->n, i);
        exit(-1);
    }

    initv(sim->h, rmax, rsub, sim->hsub, levtot, sim->rlev, sim->hlev, steps);
}
