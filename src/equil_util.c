
#include "equil_util.h"

#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "consts.h"
#include "nrutil.h"

/***************************************************************************/
/* Routine that locates nearest grid point for a given value.              */
/* Adapted from Numerical Recipes.                                         */
/* It returns jlow such that xx[jlo] ≤ x < xx[jhi]                         */
/***************************************************************************/
void hunt(double xx[], int n, double x, int *jlo) {
    int jm, jhi, inc, ascnd;
    
    ascnd = (xx[n] > xx[1]);
    if (*jlo <= 0 || *jlo > n) {
        *jlo = 0;
        jhi = n + 1;
    } else {
        inc = 1;
        if ((x >= xx[*jlo]) == ascnd) {
            if (*jlo == n) return;
            jhi = (*jlo) + 1;
            while ((x >= xx[jhi]) == ascnd) {
                *jlo = jhi;
                inc += inc;
                jhi = (*jlo) + inc;
                if (jhi > n) {
                    jhi = n + 1;
                    break;
                }
            }
        } else {
            if (*jlo == 1) {
                *jlo = 0;
                return;
            }
            jhi = (*jlo);
            *jlo -= 1;
            while ((x < xx[*jlo]) == ascnd) {
                jhi = (*jlo);
                inc += inc;
                *jlo = jhi - inc;
                if (*jlo < 1) {
                    *jlo = 0;
                    break;
                }
            }
        }
    }
    while (jhi - (*jlo) != 1) {
        jm = (jhi + (*jlo)) >> 1;
        if ((x > xx[jm]) == ascnd)
            *jlo = jm;
        else
            jhi = jm;
    }
}

/*C*/
/*************************************************************************/
/* Driver for the interpolation routine. First we find the tab. point    */
/* nearest to xb, then we interpolate around xb.                         */
/*************************************************************************/
/* Linear interpolation */
double interp(double xp[],
              double yp[],
              int np,
              double xb,
              int *n_nearest_pt) {
    int k;    /* index of 1st point */
    double y; /* intermediate value */

    hunt(xp, np, xb, n_nearest_pt);
    k = IMIN(IMAX(*n_nearest_pt, 1), np);

    if (xb == xp[k]) {
        xb += DBL_EPSILON;
    }

    y = (xb - xp[k]) * (yp[k + 1] - yp[k]) / (xp[k + 1] - xp[k]) + yp[k];

    return (y);
}

/*******************************************************************/
/* Returns the derivative w.r.t. s                                 */
/*******************************************************************/
double deriv_s(double **f, int s, int m) {
    double d_temp;

    switch (s) {
        case 1:
            d_temp = (f[s + 1][m] - f[s][m]) / DS;
            break;

        case SDIV:
            d_temp = (f[s][m] - f[s - 1][m]) / DS;
            break;

        default:
            d_temp = (f[s + 1][m] - f[s - 1][m]) / (2.0 * DS);
            break;
    }
    return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. s                                 */
/*******************************************************************/
// double deriv_ss(double **f, int s, int m) {
//     double d_temp;

//     switch (s) {
//         case 1:
//             s = 4;
//             d_temp = (f[s + 2][m] - 2.0 * f[s][m] + f[s - 2][m]) / (4.0 * SQ(DS));
//             break;

//         case 2:
//             s = 4;
//             d_temp = (f[s + 2][m] - 2.0 * f[s][m] + f[s - 2][m]) / (4.0 * SQ(DS));
//             break;

//         case 3:
//             s = 4;
//             d_temp = (f[s + 2][m] - 2.0 * f[s][m] + f[s - 2][m]) / (4.0 * SQ(DS));
//             break;

//         case SDIV - 1:
//             s = SDIV - 2;
//             d_temp = (f[s + 2][m] - 2.0 * f[s][m] + f[s - 2][m]) / (4.0 * SQ(DS));
//             break;

//         case SDIV:
//             s = SDIV - 2;
//             d_temp = (f[s + 2][m] - 2.0 * f[s][m] + f[s - 2][m]) / (4.0 * SQ(DS));
//             break;

//         default:
//             d_temp = (f[s + 2][m] - 2.0 * f[s][m] + f[s - 2][m]) / (4.0 * SQ(DS));
//             break;
//     }
//     return d_temp;
// }

double deriv_ss(double **f, int s, int m) {
    double d_temp;

    if (s >= 2 && s <= SDIV - 2) {
        // Central second-order accurate
        d_temp = (f[s + 2][m] - 2.0 * f[s][m] + f[s - 2][m]) / (4.0 * SQ(DS));
    } else if (s == 0) {
        // Forward difference
        d_temp = (2.0 * f[s][m] - 5.0 * f[s + 1][m] + 4.0 * f[s + 2][m] - f[s + 3][m]) / SQ(DS);
    } else if (s == 1) {
        // Forward-biased 2nd derivative
        d_temp = (f[s + 2][m] - 2.0 * f[s][m] + f[s - 1][m]) / SQ(DS);  // 2nd-order mixed
    } else if (s == SDIV - 1) {
        // Backward-biased
        d_temp = (f[s + 1][m] - 2.0 * f[s][m] + f[s - 2][m]) / SQ(DS);  // Mixed
    } else if (s == SDIV) {
        // Backward difference
        d_temp = (2.0 * f[s][m] - 5.0 * f[s - 1][m] + 4.0 * f[s - 2][m] - f[s - 3][m]) / SQ(DS);
    } else {
        // Safety net
        d_temp = 0.0;
    }

    return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. mu                                */
/*******************************************************************/
double deriv_m(double **f, int s, int m) {
    double d_temp;

    switch (m) {
        case 1:
            d_temp = (f[s][m + 1] - f[s][m]) / DM;
            break;

        case MDIV:
            d_temp = (f[s][m] - f[s][m - 1]) / DM;
            break;

        default:
            d_temp = (f[s][m + 1] - f[s][m - 1]) / (2.0 * DM);
            break;
    }
    return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. s                                 */
/*******************************************************************/
// double deriv_mm(double **f, int s, int m) {
//     double d_temp;

//     switch (m) {
//         case 1:
//             m = 2;
//             d_temp = (f[s][m + 1] - 2.0 * f[s][m] + f[s][m - 1]) / SQ(DM);
//             break;

//         case MDIV:
//             m = MDIV - 1;
//             d_temp = (f[s][m + 1] - 2.0 * f[s][m] + f[s][m - 1]) / SQ(DM);
//             break;

//         default:
//             d_temp = (f[s][m + 1] - 2.0 * f[s][m] + f[s][m - 1]) / SQ(DM);
//             break;
//     }
//     return d_temp;
// }

double deriv_mm(double **f, int s, int m) {
    double d_temp;

    if (m >= 2 && m <= MDIV - 1) {
        // Central second-order
        d_temp = (f[s][m + 1] - 2.0 * f[s][m] + f[s][m - 1]) / SQ(DM);
    } else if (m == 1) {
        // Forward-biased second derivative
        d_temp = (2.0 * f[s][m] - 5.0 * f[s][m + 1] + 4.0 * f[s][m + 2] - f[s][m + 3]) / SQ(DM);
    } else if (m == MDIV) {
        // Backward-biased second derivative
        d_temp = (2.0 * f[s][m] - 5.0 * f[s][m - 1] + 4.0 * f[s][m - 2] - f[s][m - 3]) / SQ(DM);
    } else {
        // Failsafe: outside valid domain
        d_temp = 0.0;
    }
    return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. s and mu                          */
/*******************************************************************/
double deriv_sm(double **f, int s, int m) {
    double d_temp;

    switch (s) {
        case 1:
            if (m == 1) {
                d_temp = (f[s + 1][m + 1] - f[s][m + 1] - f[s + 1][m] + f[s][m]) / (DM * DS);
            } else {
                if (m == MDIV) {
                    d_temp = (f[s + 1][m] - f[s][m] - f[s + 1][m - 1] + f[s][m - 1]) / (DM * DS);
                } else {
                    d_temp = (f[s + 1][m + 1] - f[s + 1][m - 1] - f[s][m + 1] + f[s][m - 1]) /
                             (2.0 * DM * DS);
                }
            }
            break;

        case SDIV:
            if (m == 1) {
                d_temp = (f[s][m + 1] - f[s][m] - f[s - 1][m + 1] + f[s - 1][m]) / (DM * DS);
            } else {
                if (m == MDIV) {
                    d_temp = (f[s][m] - f[s - 1][m] - f[s][m - 1] + f[s - 1][m - 1]) / (DM * DS);
                } else {
                    d_temp = (f[s][m + 1] - f[s][m - 1] - f[s - 1][m + 1] + f[s - 1][m - 1]) /
                             (2.0 * DM * DS);
                }
            }
            break;

        default:
            if (m == 1) {
                d_temp = (f[s + 1][m + 1] - f[s - 1][m + 1] - f[s + 1][m] + f[s - 1][m]) / (2.0 * DM * DS);
            } else {
                if (m == MDIV) {
                    d_temp = (f[s + 1][m] - f[s - 1][m] - f[s + 1][m - 1] + f[s - 1][m - 1]) /
                             (2.0 * DM * DS);
                } else {
                    d_temp = (f[s + 1][m + 1] - f[s - 1][m + 1] - f[s + 1][m - 1] + f[s - 1][m - 1]) /
                             (4.0 * DM * DS);
                }
            }
            break;
    }

    return d_temp;
}

/*******************************************************************/
/* Returns the Legendre polynomial of degree n, evaluated at x.    */
/*******************************************************************/
double legendre(int n, double x) /* checked */
{
    int i; /* counter */

    double p, /* Legendre polynomial of order n */
        p_1,  /*    "         "      "    "   n-1*/
        p_2;  /*    "         "      "    "   n-2 */

    p_2 = 1.0;
    p_1 = x;

    if (n >= 2) {
        for (i = 2; i <= n; i++) {
            p = (x * (2.0 * i - 1.0) * p_1 - (i - 1.0) * p_2) / i;
            p_2 = p_1;
            p_1 = p;
        }
        return p;
    } else {
        if (n == 1)
            return p_1;
        else
            return p_2;
    }
}

/*******************************************************************/
/* Returns the associated Legendre polynomial P_l^m(x).            */
/* Adapted from numerical recipes.                                 */
/*******************************************************************/
double plgndr(int l, int m, double x) {
    double fact, pll, pmm, pmmp1, somx2;
    int i, ll;

    if (m < 0 || m > l || fabs(x) > 1.0)
        printf("Bad arguments in routine PLGNDR");
    pmm = 1.0;
    if (m > 0) {
        somx2 = sqrt((1.0 - x) * (1.0 + x));
        fact = 1.0;
        for (i = 1; i <= m; i++) {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    if (l == m)
        return pmm;
    else {
        pmmp1 = x * (2 * m + 1) * pmm;
        if (l == (m + 1))
            return pmmp1;
        else {
            for (ll = (m + 2); ll <= l; ll++) {
                pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
                pmm = pmmp1;
                pmmp1 = pll;
            }
            return pll;
        }
    }
}

void compute_polynomials(double *s_gp, double *mu,
                         float ***f_rho, float ***f_gama,
                         double **P_2n, double **P1_2n_1,
                         double *sin_theta, double *theta) {
    int n, i, j, k, m;
    volatile double sk, sj, sk1, sj1;

    double **f2n = dmatrix(1, LMAX + 1, 1, SDIV);
    /* Compute f2n */
    for (n = 0; n <= LMAX; n++) {
        for (i = 2; i <= SDIV; i++) {
            f2n[n + 1][i] = pow((1.0 - s_gp[i]) / s_gp[i], 2.0 * n);
        }
    }

    /* Compute f_rho and f_gama */
    if (SMAX != 1.0) {
        for (j = 2; j <= SDIV; j++) {
            for (n = 1; n <= LMAX; n++) {
                for (k = 2; k <= SDIV; k++) {
                    sk = s_gp[k];
                    sj = s_gp[j];
                    sk1 = 1.0 - sk;
                    sj1 = 1.0 - sj;
                    if (k < j) {
                        f_rho[j][n + 1][k] = f2n[n + 1][j] * sj1 /
                                             (sj * f2n[n + 1][k] * sk1 * sk1);
                        f_gama[j][n + 1][k] = f2n[n + 1][j] /
                                              (f2n[n + 1][k] * sk * sk1);
                    } else {
                        f_rho[j][n + 1][k] = f2n[n + 1][k] /
                                             (f2n[n + 1][j] * sk * sk1);
                        f_gama[j][n + 1][k] = f2n[n + 1][k] * sj1 * sj1 * sk /
                                              (sj * sj * f2n[n + 1][j] * sk1 * sk1 * sk1);
                    }
                }
            }
        }
        j = 1;
        n = 0;
        for (k = 2; k <= SDIV; k++) {
            sk = s_gp[k];
            f_rho[j][n + 1][k] = 1.0 / (sk * (1.0 - sk));
        }
        n = 1;
        for (k = 2; k <= SDIV; k++) {
            sk = s_gp[k];
            sk1 = 1.0 - sk;
            f_rho[j][n + 1][k] = 0.0;
            f_gama[j][n + 1][k] = 1.0 / (sk * sk1);
        }
        for (n = 2; n <= LMAX; n++) {
            for (k = 1; k <= SDIV; k++) {
                f_rho[j][n + 1][k] = 0.0;
                f_gama[j][n + 1][k] = 0.0;
            }
        }
        k = 1;
        n = 0;
        for (j = 1; j <= SDIV; j++)
            f_rho[j][n + 1][k] = 0.0;
        for (j = 1; j <= SDIV; j++) {
            for (n = 1; n <= LMAX; n++) {
                f_rho[j][n + 1][k] = 0.0;
                f_gama[j][n + 1][k] = 0.0;
            }
        }
        n = 0;
        for (j = 2; j <= SDIV; j++) {
            for (k = 2; k <= SDIV; k++) {
                sk = s_gp[k];
                sj = s_gp[j];
                sk1 = 1.0 - sk;
                sj1 = 1.0 - sj;
                if (k < j)
                    f_rho[j][n + 1][k] = sj1 / (sj * sk1 * sk1);
                else
                    f_rho[j][n + 1][k] = 1.0 / (sk * sk1);
            }
        }
    } else {
        for (j = 2; j <= SDIV - 1; j++) {
            for (n = 1; n <= LMAX; n++) {
                for (k = 2; k <= SDIV - 1; k++) {
                    sk = s_gp[k];
                    sj = s_gp[j];
                    sk1 = 1.0 - sk;
                    sj1 = 1.0 - sj;
                    if (k < j) {
                        f_rho[j][n + 1][k] = f2n[n + 1][j] * sj1 /
                                             (sj * f2n[n + 1][k] * sk1 * sk1);
                        f_gama[j][n + 1][k] = f2n[n + 1][j] /
                                              (f2n[n + 1][k] * sk * sk1);
                    } else {
                        f_rho[j][n + 1][k] = f2n[n + 1][k] /
                                             (f2n[n + 1][j] * sk * sk1);
                        f_gama[j][n + 1][k] = f2n[n + 1][k] * sj1 * sj1 * sk /
                                              (sj * sj * f2n[n + 1][j] * sk1 * sk1 * sk1);
                    }
                }
            }
        }
        j = 1;
        n = 0;
        for (k = 2; k <= SDIV - 1; k++) {
            sk = s_gp[k];
            f_rho[j][n + 1][k] = 1.0 / (sk * (1.0 - sk));
        }
        n = 1;
        for (k = 2; k <= SDIV - 1; k++) {
            sk = s_gp[k];
            sk1 = 1.0 - sk;
            f_rho[j][n + 1][k] = 0.0;
            f_gama[j][n + 1][k] = 1.0 / (sk * sk1);
        }
        for (n = 2; n <= LMAX; n++) {
            for (k = 1; k <= SDIV - 1; k++) {
                f_rho[j][n + 1][k] = 0.0;
                f_gama[j][n + 1][k] = 0.0;
            }
        }
        k = 1;
        n = 0;
        for (j = 1; j <= SDIV - 1; j++)
            f_rho[j][n + 1][k] = 0.0;
        for (j = 1; j <= SDIV - 1; j++) {
            for (n = 1; n <= LMAX; n++) {
                f_rho[j][n + 1][k] = 0.0;
                f_gama[j][n + 1][k] = 0.0;
            }
        }
        n = 0;
        for (j = 2; j <= SDIV - 1; j++) {
            for (k = 2; k <= SDIV - 1; k++) {
                sk = s_gp[k];
                sj = s_gp[j];
                sk1 = 1.0 - sk;
                sj1 = 1.0 - sj;
                if (k < j)
                    f_rho[j][n + 1][k] = sj1 / (sj * sk1 * sk1);
                else
                    f_rho[j][n + 1][k] = 1.0 / (sk * sk1);
            }
        }
        j = SDIV;
        for (n = 1; n <= LMAX; n++) {
            for (k = 1; k <= SDIV; k++) {
                f_rho[j][n + 1][k] = 0.0;
                f_gama[j][n + 1][k] = 0.0;
            }
        }
        k = SDIV;
        for (j = 1; j <= SDIV; j++) {
            for (n = 1; n <= LMAX; n++) {
                f_rho[j][n + 1][k] = 0.0;
                f_gama[j][n + 1][k] = 0.0;
            }
        }
    }

    /* Compute Legendre functions */
    n = 0;
    for (i = 1; i <= MDIV; i++) {
        P_2n[i][n + 1] = legendre(2 * n, mu[i]);
    }
    for (i = 1; i <= MDIV; i++) {
        for (n = 1; n <= LMAX; n++) {
            P_2n[i][n + 1] = legendre(2 * n, mu[i]);
            P1_2n_1[i][n + 1] = plgndr(2 * n - 1, 1, mu[i]);
        }
    }

    /* Free temporary matrix */
    free_dmatrix(f2n, 1, LMAX + 1, 1, SDIV);

    /* Allocate and compute sin_theta and theta arrays */
    for (m = 1; m <= MDIV; m++) {
        sin_theta[m] = sqrt(1.0 - mu[m] * mu[m]);
        theta[m] = asin(sin_theta[m]);
    }
}

/*******************************************************************/
// ROTATION LAWS
/* J constant law */

double call_diff_rotation(double first_arg, double *params) {
    return diff_rotation(first_arg, params[0], params[1], params[2], params[3],
                         params[4], params[5], params[6], params[7]);
}
double diff_rotation(double x,
                     double r_e,
                     double s_e,
                     double rho_equator_h,
                     double gama_equator_h,
                     double omega_equator_h,
                     double rho_pole_h,
                     double gama_pole_h,
                     double A_diff) {
    double v = (x - omega_equator_h) * (s_e / (1 - s_e)) * exp(-SQ(r_e) * rho_equator_h);
    return (SQ(r_e) * (gama_equator_h + rho_equator_h - gama_pole_h - rho_pole_h) + log(1.0 - SQ(v))) * SQ(1.0 - SQ(v)) - SQ((1.0 / A_diff) * (x - omega_equator_h) * SQ(s_e / (1 - s_e)) * exp(-2.0 * SQ(r_e) * rho_equator_h));
}

double call_rotation_law(double first_arg, double *params) {
    return rotation_law(first_arg, params[0], params[1], params[2], params[3], params[4], params[5],
                        params[6]);
}
double rotation_law(double x,
                    double r_e,
                    double rhogp,
                    double omegagp,
                    double sgp,
                    double mugp,
                    double Omega_c,
                    double A_diff) {
    return (Omega_c - x) * (SQ(1.0 - sgp) - SQ(sgp * (x - omegagp) * exp(-SQ(r_e) * rhogp)) * (1.0 - mugp * mugp)) - SQ((1.0 / A_diff)) * (x - omegagp) * sgp * sgp * (1.0 - mugp * mugp) * exp(-2.0 * SQ(r_e) * rhogp);
}

/* Uryu 8, Eq. (16) of https://arxiv.org/pdf/2405.06609 */

double call_uryu_diff_rotation(double first_arg, double *params) {
    return uryu_diff_rotation(first_arg, params[0], params[1], params[2], params[3], params[4], params[5],
                              params[6], params[7], params[8], params[9]);
}
double uryu_diff_rotation(double x,
                          double r_e,
                          double s_e,
                          double rho_equator_h,
                          double gama_equator_h,
                          double omega_equator_h,
                          double rho_pole_h,
                          double gama_pole_h,
                          double A_diff, double B_diff,
                          double lambda2) {
    /* EQ A7 of  https://doi.org/10.1093/mnras/stab392*/

    // Auxiliary variables
    double Omo = x - omega_equator_h;  // \Omega_eq_hat - \omega_eq_hat
    double exprho2 = exp(-2. * SQ(r_e) * rho_equator_h);
    double D = 1. - SQ(Omo) * SQ(s_e / (1. - s_e)) * exprho2;
    double Je = Omo * SQ(s_e) * exprho2;
    Je /= (SQ(1. - s_e) - Omo * Je);

    double p1 = SQ(r_e) * (gama_equator_h + rho_equator_h - gama_pole_h - rho_pole_h) + log(D);
    double p2 = 2. * x * Je;
    double p3 = 2. * SQ(A_diff / B_diff) * atan2(SQ(lambda2 * Je), SQ(SQ(A_diff) * x));
    p3 -= sqrt(2.) * (atan2(SQ(A_diff) * x - lambda2 * sqrt(2.) * Je, SQ(A_diff) * x) - atan2(SQ(A_diff) * x + lambda2 * sqrt(2.) * Je, SQ(A_diff) * x));

    double N1 = (SQ(A_diff) * x * sqrt(2.) / lambda2) * Je;
    double N2 = SQ(Je) + SQ(SQ(A_diff) * x / lambda2);
    p3 += sqrt(2.) * atanh(N1 / N2);

    return p1 + p2 - SQ(A_diff * x / lambda2) / 2.   * p3;
}

double call_uryu_rotation_law(double first_arg, double *params) {
    return uryu_rotation_law(first_arg, params[0], params[1], params[2], params[3], params[4], params[5],
                             params[6], params[7]);
}
double uryu_rotation_law(double x,
                         double r_e,
                         double rhogp,
                         double omegagp,
                         double sgp,
                         double mugp,
                         double Omega_c,
                         double A_diff, double B_diff) {
    /* EQ A8 of  https://doi.org/10.1093/mnras/stab392 for p = 1, q = 3*/

    double const p = 1.;
    double const q = 3.;

    double J = (x - omegagp) * SQ(sgp) * (1 - SQ(mugp)) * exp(-2. * SQ(r_e) * rhogp);
    J /= (SQ(1. - sgp) - (x - omegagp) * J);
    double N = 1. + pow(J / SQ(B_diff) / Omega_c, p);
    double D = 1. + pow(J / SQ(A_diff) / Omega_c, p + q);

    return x * D - Omega_c * N;
}

// Compute J_max
double call_uryu_dOdj(double first_arg, double *params) {
    return uryu_dOdj(first_arg, params[0], params[1], params[2], params[3], params[4]);
}
double uryu_dOdj(double j, double Oc, double p, double q, double A, double B) {
    double x1 = pow(j / (SQ(A) * Oc), p + q);
    double x2 = pow(j / (SQ(B) * Oc), p);

    double num = p * (1 + x1) * pow(x2, 1 - 1 / p) / SQ(B) - (p + q) * (1 + x2) * pow(x1, 1 - 1 / (p + q)) / SQ(A);
    double den = SQ(1 + x1);
    return num / den;
}

/* Uryu Extended, Eq. (31) of https://arxiv.org/pdf/2405.06609 */
// Compute J_max
double call_uryuExt_dOdj(double first_arg, double *params) {
    return uryuExt_dOdj(first_arg, params[0], params[1], params[2], params[3], params[4]);
}
double uryuExt_dOdj(double j, double Oc, double p, double csi, double A, double B) {
    return (p * pow(j / Oc, p) - pow(A, -2 * csi) * (pow(B, 2 * p) * csi + (p + csi) * pow(j / Oc, p)) * pow(j / Oc, csi)) / j;
}

// Compute B^2p Omega_c^p
double call_uryuExt_Beq(double first_arg, double *params) {
    return uryuExt_Beq(first_arg, params[0], params[1], params[2], params[3], params[4], params[5]);
}
double uryuExt_Beq(double x, double p, double csi, double l1, double l2, double je, double jm) {
    /* x = B^2p Omega_c^p*/
    double alp = pow(jm, csi) * (1 - l2) - pow(je, csi) * (1 - l1);
    double bet = pow(jm, p + csi) * (1 - l2) - pow(je, csi) * pow(jm, p) + pow(je, p) * pow(jm, csi) - pow(je, p + csi) * (1 - l1);
    double gam = pow(jm, p + csi) * pow(je, p) - pow(je, p + csi) * pow(jm, p);

    return alp * x * x + bet * x + gam;
}

// Compute Omega_e
double call_uryuExt_diff_rotation(double first_arg, double *params) {
    return uryuExt_diff_rotation(first_arg, params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8], params[9], params[10], params[11]);
}
double uryuExt_diff_rotation(double x,
                             double s_e,
                             double r_e,
                             double rho_equator_h,
                             double gama_equator_h,
                             double omega_equator_h,
                             double rho_pole_h,
                             double gama_pole_h,
                             double A_diff, double B_diff,
                             double lambda2, double p9, double csi) {
    /* EQ 39 + 50 of https://arxiv.org/pdf/2405.06609 */

    // Auxiliary variables
    double Oc = x / lambda2;
    double Omo = x - omega_equator_h;  // \Omega_eq_hat - \omega_eq_hat
    double Omo2 = SQ(Omo);
    double re2 = SQ(r_e);
    double expnu2 = exp(-2. * re2 * rho_equator_h);
    double D = 1. - Omo2 * SQ(s_e) / SQ(1 - s_e) * expnu2;

    double Je = Omo * SQ(s_e) * expnu2 / (SQ(1 - s_e) - Omo2 * SQ(s_e) * expnu2);

    double p1 = re2 * (gama_equator_h + rho_equator_h - gama_pole_h - rho_pole_h) + log(D);
    double p2 = uryuExt_integral(Je, Oc, p9, csi, A_diff, B_diff);

    return p1 + 2 * Oc * p2;
}

// Compute \int j \partial{Omega} / \partial{j} dj
double uryuExt_integral(double J, double Oc, double p, double csi, double A, double B) {
    double p2 = (p / (p + 1)) * pow(J, p + 1) / pow(SQ(B) * Oc, p);
    p2 += -(p / (p + csi + 1)) * pow(J, p + csi + 1) / (pow(SQ(B) * Oc, p) * pow(SQ(A) * Oc, csi));
    p2 += -(csi / (csi + 1)) * pow(J, csi + 1) / pow(SQ(A) * Oc, csi);
    p2 += -(csi / (p + csi + 1)) * pow(J, p + csi + 1) / (pow(SQ(B) * Oc, p) * pow(SQ(A) * Oc, csi));

    return p2;
}

// Compute Omega(j)
double call_uryuExt_rotation_law(double first_arg, double *params) {
    return uryuExt_rotation_law(first_arg, params[0], params[1], params[2], params[3], params[4], params[5],
                                params[6], params[7], params[8], params[9]);
}
double uryuExt_rotation_law(double x,
                            double r_e,
                            double rhogp,
                            double omegagp,
                            double sgp,
                            double mugp,
                            double Omega_c,
                            double A_diff, double B_diff,
                            double p9, double csi) {
    double n1 = (x - omegagp) * SQ(sgp) * (1 - SQ(mugp)) * exp(-2. * SQ(r_e) * rhogp);
    double n2 = (SQ(1. - sgp) - SQ(x - omegagp) * SQ(sgp) * (1. - SQ(mugp)) * exp(-2. * SQ(r_e) * rhogp));
    double j = n1 / n2;

    double tmp = (1 + pow(j / (SQ(B_diff) * Omega_c), p9)) * (1 - pow(j / (SQ(A_diff) * Omega_c), csi));

    return x - Omega_c * tmp;
}

////////////////
#include <stdint.h>
int my_isnan(double x) {
    uint64_t bits;
    memcpy(&bits, &x, sizeof(x));
    uint64_t exp  = (bits >> 52) & 0x7FF;
    uint64_t frac = bits & 0xFFFFFFFFFFFFF;
    return (exp == 0x7FF && frac != 0);
}
// maximum = true is used find the value x0 such that func(x = x0) := dF/dx|_x0 = 0 and F(x0) is a local maximum of the function (i.e. it is used to find the Jm_h such that dOmega/dj|_Jm_h = 0)
int find_boundary(double (*func)(double, double *), double start, double step, double *inf, double *sup, bool maximum, bool stop, char *varName, int num_args, ...) {
    va_list args;
    va_start(args, num_args);

    double *arg_array;

    arg_array = (double *)malloc(num_args * sizeof(double));
    if (!arg_array) {
        printf("Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < num_args; i++) {
        arg_array[i] = va_arg(args, double);
    }
    va_end(args);

    double xp = 0., yp = 0., x = 0., y = 0.;
    int iter = 0;
    int fail = 0;

    x = start;
    y = (*func)(x, arg_array);
    double micro_step = step < 0 ? - 0.001 : 0.001;

    // if (strcmp(varName, "Omega_e") == 0) 
    //     printf("%g  %g , %g   %g, %g\n", xp, yp, x, y, micro_step);

    while (my_isnan(y) && iter < 1000) {
        x += micro_step;
        y = (*func)(x, arg_array);
        iter++;
    };

    // if (strcmp(varName, "Omega_e") == 0) 
    //     printf("%g  %g , %g   %g, %g\n", xp, yp, x, y, micro_step);

    double x0 = x;
    double y0 = y;

    for (int i = 0; i <= NTRIES; i++) {
        double tmp_step = step * pow(10, -i);
        x = x0;
        y = y0;

        iter = 0;
        fail = 1;
        while (iter < pow(10., (double)i) * MAXIT) {
            xp = x;
            yp = y;
            x += tmp_step;
            y = (*func)(x, arg_array);

            if (my_isnan(y) && my_isnan(yp)) {
                break;
            }

            if (y * yp <= 0 && (!maximum || (yp > 0 && y < 0))) {
                fail = 0;
                break;
            }

            // if (strcmp(varName, "Omega_e") == 0) 
            //     printf("%g  %g , %g   %g\n", xp, yp, x, y);

            iter++;
        }
        if (fail == 0) {
            break;
        }
    }
    // if (strcmp(varName, "Omega_e") == 0) 
    //             printf("%g  %g , %g   %g\n", xp, yp, x, y);
    if (fail == 1) {
        if (stop) {
            printf("Failed to find brackets for %s!\n", varName);
            exit(EXIT_FAILURE);
        } else {
            return 1;
        }
    }

    *inf = xp;
    *sup = x;
    free(arg_array);
    return 0;
}

/* From `Numerical Recipies' chapter 9.3 */
double brent_root_finder(double (*func)(double, double *),
                         double a, double b,
                         double tol, char *varName, int num_args, ...) {
    int iter;
    double c,
        d,
        e,
        min1,
        min2;

    double fa, fb, fc,
        p,
        q,
        r,
        s,
        tol1,
        xm;

    va_list args;
    va_start(args, num_args);

    double *arg_array;

    arg_array = (double *)malloc(num_args * sizeof(double));
    if (!arg_array) {
        printf("Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < num_args; i++) {
        arg_array[i] = va_arg(args, double);
    }
    va_end(args);

    fa = (*func)(a, arg_array);
    fb = (*func)(b, arg_array);

    double x1 = a;
    double x2 = b;
    double fx1 = fa;
    double fx2 = fb;

    c = b;
    fc = fb;
    for (iter = 1; iter <= MAXIT; iter++) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a;
            fc = fa;
            e = d = b - a;
        }

        if (fabs(fc) < fabs(fb)) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol;
        xm = 0.5 * (c - b);

        if (fabs(xm) <= tol1 || fabs(fb) <= EPS) {
            free(arg_array);
            return b;
        }

        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            s = fb / fa;
            if (a == c) {
                // Linear interpolation (special case)
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                // Quadratic interpolation formula
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) {
                q = -q;
            }
            p = fabs(p);
            min1 = 3.0 * xm * q - fabs(tol1 * q);
            min2 = fabs(e * q);

            if (2.0 * p < fmin(min1, min2)) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }

        a = b;
        fa = fb;

        if (fabs(d) > tol1) {
            b += d;
        } else {
            b += SIGN(tol1, xm);
        }
        fb = (*func)(b, arg_array);
    }

    printf("Maximum number of iterations exceeded in ZBRENT for %s!\n", varName);
    printf("Upper bound %.4e at x %.4e\n", fx2, x2);
    printf("Lower bound %.4e at x %.4e\n", fx1, x1);
    fflush(stdout);

    free(arg_array);
    exit(EXIT_FAILURE);
}