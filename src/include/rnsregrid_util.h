#include "consts.h"
#include "hdf5.h"

#ifndef _RNSREGRID_UTILS_H_
#define _RNSREGRID_UTILS_H_

double **array_allocate(long nrl, long nrh, long ncl, long nch);
void array_free(double **m, long nrl, long nrh, long ncl, long nch);

double ***tensor_allocate(long nrl,
                          long nrh,
                          long ncl,
                          long nch,
                          long ndl,
                          long ndh);

void tensor_free(double ***t,
                 long nrl,
                 long nrh,
                 long ncl,
                 long nch,
                 long ndl,
                 long ndh);

void calculate_chunk_dims(hsize_t *chunk_dims, hsize_t *dims);

double interp_4(double xp[5],
                double yp[5],
                int np,
                double xb);

void grid_interp(double **old,
                 double *s_gp,
                 double *mu,
                 double r_e,
                 int nx,
                 int ny,
                 int nz,
                 double *x_grid,
                 double *y_grid,
                 double *z_grid,
                 int i,
                 int j,
                 int k,
                 double *new,
                 int sign);

void grid_interp_all(double *s_gp,
                     double *mu,
                     double r_e,
                     int nx,
                     int ny,
                     int nz,
                     double *x_grid,
                     double *y_grid,
                     double *z_grid,
                     int i,
                     int j,
                     int k,
                     double **nu,
                     double **B,
                     double **alpha,
                     double **omega,
                     double **nu_dr,
                     double **B_dr,
                     double **alpha_dr,
                     double **omega_dr,
                     double **nu_dth,
                     double **B_dth,
                     double **alpha_dth,
                     double **omega_dth,
                     double **rho_0_BM,
                     double **energy_BM,
                     double **pressure_BM,
                     double **rho_0_DM,
                     double **energy_DM,
                     double **pressure_DM,
                     double **Omega_BM,
                     double **Omega_DM,
                     double *nu_c,
                     double *B_c,
                     double *alpha_c,
                     double *omega_c,
                     double *nu_dr_c,
                     double *B_dr_c,
                     double *alpha_dr_c,
                     double *omega_dr_c,
                     double *nu_dth_c,
                     double *B_dth_c,
                     double *alpha_dth_c,
                     double *omega_dth_c,
                     double *rho_0_BM_c,
                     double *energy_BM_c,
                     double *pressure_BM_c,
                     double *rho_0_DM_c,
                     double *energy_DM_c,
                     double *pressure_DM_c,
                     double *distance_c,
                     double *OmegaBM_c,
                     double *OmegaDM_c);

#endif
