#ifndef MAIN_UTIL_H
#define MAIN_UTIL_H

#include <stdbool.h>
#include "consts.h"

void set_default_values(double *cf, double *accuracy, double *r_ratio_BM, double *r_ratio_DM, double *fdm_target, double *MDM, double *fdm_tol,  int *verbose, struct stellar_properties *star_props, struct stellar_properties *DM_props, struct EOS *eosBM, struct EOS *eosDM, struct DiffRotParams *DiffRotBM, struct DiffRotParams *DiffRotDM, char *input_filename, char *output2d_filename, char *outputS_filename, char *output_name, char *id_file);

void outfiles_name(char *output2d_filename, char *outputS_filename, char *output_name, size_t bufsize);

void zero_DM(struct evolution_variables *ev, struct stellar_properties *DM_props, double *r_ratio_DM_target, double *r_e_DM);

void copy_fields(struct evolution_variables *Dest, struct evolution_variables *Src);

void load_initial_data(char *id_file, double *s_gp, double *mu, struct evolution_variables *evolution_functions, struct evolution_variables *tmp_func, struct stellar_properties *star_props, struct stellar_properties *DM_props, struct EOS *eosBM, struct EOS *eosDM, double *r_ratio_BM, double *r_ratio_DM, double *Omega_BM, double *Omega_DM, double *r_e_BM, double *r_e_DM, double *s_e_BM, double *s_e_DM, struct DiffRotParams *DiffRotBM, struct DiffRotParams *DiffRotDM, bool counter, double *Mass_BM, double *Mass_0_BM, double *J_BM, double *R_e_BM, double *v_plus, double *v_minus, double *Omega_K_BM, double *Mass_DM, double *Mass_0_DM, double *J_DM,  double *R_e_DM, double MDM, double *R_p_BM, double *R_p_DM, double *Omega_K_DM, double *T_BM, double *T_DM, double *W_BM, double *W_DM, double *GRV2, double *GRV3, int verbose);

void regrid(double *s_gp, double *tmp_sgp, int tmp_sdiv, double *mu, double *tmp_mu, int tmp_mdiv, struct evolution_variables *ev, struct evolution_variables *tmp_ev);

void do_jconstant(double *s_gp, double *mu, struct evolution_variables *ev, struct EOS *eosBM, struct EOS *eosDM, struct stellar_properties *star_props, struct stellar_properties *DM_props, struct DiffRotParams *DiffRotBM, struct DiffRotParams *DiffRotDM, double old_rrBM, double old_rrDM, double *r_e_BM, double *Omega_BM, double *s_e_BM, double *r_e_DM, double *Omega_DM, double *s_e_DM, bool *outOfiter, double *Omega_e_BM, double *Omega_e_DM, double r_ratio_BM, double r_ratio_DM, bool *BM_first_diff_done, bool *DM_first_diff_done, bool zero, double accuracy, double cf, bool counter, int verbose);

#endif