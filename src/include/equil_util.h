#include <stdbool.h>

#ifndef EQUIL_UTIL_H
#define EQUIL_UTIL_H

void hunt(double xx[], int n, double x, int *jlo);  // EDO: looks for neighb. points

double interp(double xp[],
              double yp[],
              int np,
              double xb,
              int *n_nearest_pt);  // interpolation //it calls hunt

double deriv_s(double **f, int s, int m);  // all kind of derivatives

double deriv_ss(double **f, int s, int m);

double deriv_m(double **f, int s, int m);

double deriv_mm(double **f, int s, int m);

double deriv_sm(double **f, int s, int m);

double legendre(int n, double x);

double plgndr(int l, int m, double x);

void compute_polynomials(double *s_gp, double *mu,
                         float ***f_rho, float ***f_gama,
                         double **P_2n, double **P1_2n_1,
                         double *sin_theta, double *theta);

/*******************************************************************/

/* J constant */
double call_diff_rotation(double first_arg, double *params);
double diff_rotation(double x,
                     double r_e,
                     double s_e,
                     double rho_equator_h,
                     double gama_equator_h,
                     double omega_equator_h,
                     double rho_pole_h,
                     double gama_pole_h,
                     double A_diff);

double call_rotation_law(double first_arg, double *params);
double rotation_law(double x,
                    double r_e,
                    double rhogp,
                    double omegagp,
                    double sgp,
                    double mugp,
                    double Omega_c,
                    double A_diff);

/* Uryu 8 */
double call_uryu_diff_rotation(double first_arg, double *params);
double uryu_diff_rotation(double x,
                          double r_e,
                          double s_e,
                          double rho_equator_h,
                          double gama_equator_h,
                          double omega_equator_h,
                          double rho_pole_h,
                          double gama_pole_h,
                          double A_diff,
                          double B_diff,
                          double lambda2);

double call_uryu_rotation_law(double first_arg, double *params);
double uryu_rotation_law(double x,
                         double r_e,
                         double rhogp,
                         double omegagp,
                         double sgp,
                         double mugp,
                         double Omega_c,
                         double A_diff,
                         double B_diff);

double call_uryu_dOdj(double first_arg, double *params);
double uryu_dOdj(double j, double Oc, double p, double q, double A, double B);

/* Uryu Extended */
double call_uryuExt_dOdj(double first_arg, double *params);
double uryuExt_dOdj(double j, double Oc, double p, double csi, double A, double B);
double call_uryuExt_Beq(double first_arg, double *params);
double uryuExt_Beq(double x, double p, double csi, double l1, double l2, double je, double jm);
double call_uryuExt_diff_rotation(double first_arg, double *params);
double uryuExt_diff_rotation(double x,
                             double s_e,
                             double r_e,
                             double rho_equator_h,
                             double gama_equator_h,
                             double omega_equator_h,
                             double rho_pole_h,
                             double gama_pole_h,
                             double A_diff, double B_diff,
                             double lambda2, double p9, double csi);
double call_uryuExt_rotation_law(double first_arg, double *params);
double uryuExt_rotation_law(double x,
                            double r_e,
                            double rhogp,
                            double omegagp,
                            double sgp,
                            double mugp,
                            double Omega_c,
                            double A_diff, double B_diff,
                            double p9, double csi);
double uryuExt_integral(double J, double Oc, double p, double csi, double A, double B);

/* Root finding */

int find_boundary(double (*func)(double, double *),
                  double start, double step, double *inf, double *sup,
                  bool maximum, bool stop, char *varName, int num_args, ...);
double brent_root_finder(double (*func)(double, double *),
                         double a, double b,
                         double tol, char *varName, int num_args, ...);

#endif