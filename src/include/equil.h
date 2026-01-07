#ifndef EQUIL_H
#define EQUIL_H

#include <stdbool.h>
#include "consts.h"

void allocate_evolution_variables(struct evolution_variables *ev);
void deallocate_evolution_variables(struct evolution_variables *ev);

void print_help();

void make_grid(double s_gp[SDIV + 1], double mu[MDIV + 1]);

void load_eos(struct EOS *eos);

double e_of_rho0(double rho0, double Gamma_P);

double e_at_p(double pp, struct EOS *eos, int *n_nearest_pt);

double p_at_e(double ee, struct EOS *eos, int *n_nearest_pt);

double p_at_h(double hh, struct EOS *eos, int *n_nearest_pt);

double h_at_p(double pp, struct EOS *eos, int *n_nearest_pt);

double n0_at_e(double ee, struct EOS *eos, int *n_nearest_pt);

double e_of_rho0_wrapper(double rho0, double *args);

void make_center(struct EOS *eos, double e_center, double *p_center, double *h_center);

void mass_radius(double *s_gp, double *mu, struct EOS *eosBM, struct evolution_variables *ev, double r_ratio_BM, double r_e_BM, double Omega_BM, double *Mass_BM, double *Mass_0_BM, double *ang_mom_BM, double *R_e_BM, double *v_plus, double *v_minus, double *Omega_K_BM, struct stellar_properties *star_props, struct EOS *eosDM, double r_ratio_DM, double r_e_DM, double Omega_DM, double *Mass_DM, double *Mass_0_DM, double *ang_mom_DM, double *R_e_DM, struct stellar_properties *DM_props, double MDM, double *R_p_BM, double *R_p_DM, double *Omega_K_DM, double *T_BM, double *T_DM, double *W_BM, double *W_DM, bool counter, double *GRV2, double *GRV3);

void compute_virial(struct evolution_variables *ev, double *s_gp, double *mu, double r_e, double r_ratio, double *grv2, double *grv3);

double dm_dr_is(double r_is, double r, double m, double p, struct stellar_properties *star_props, struct EOS *eos, int *n_nearest_pt);

double dp_dr_is(double r_is, double r, double m, double p, double ptot, struct stellar_properties *star_props, struct stellar_properties *DM_props, struct EOS *eos, int *n_nearest_pt);

double dr_dr_is(double r_is, double r, double m);

void rk4_step(double hdiv, double r_is, double r, double m, double p, double m_DM, double p_DM, bool skip_BM, bool skip_DM, struct stellar_properties *star_props, struct stellar_properties *DM_props, struct EOS *eosBM, struct EOS *eosDM, struct rk_coeff *rkc_old, struct rk_coeff *rkc_new);

void TOV(int i_check, struct stellar_properties *star_props, struct EOS *eosBM, double r_is_gp[RDIV + 1], double lambda_gp[RDIV + 1], double nu_gp[RDIV + 1], double *r_is_final_big, double *r_final_big, double *m_final, struct stellar_properties *DM_props, struct EOS *eosDM, double *r_is_final_small, double *r_final_small, double *m_DM_final, bool *is_DM);

void sphere(double s_gp[SDIV + 1], struct EOS *eosBM, struct stellar_properties *star_props, struct evolution_variables *ev, double *r_e_BM, double *s_e_BM, struct EOS *eosDM, struct stellar_properties *DM_props, double *r_e_DM, double *s_e_DM, double *massBM, double *massDM);

void computeOmega(double s_gp[SDIV + 1], double mu[MDIV + 1], double **OmegaDiff, struct evolution_variables *ev, struct stellar_properties *props, struct DiffRotParams *DiffRotParam, double r_ratio, double *omega_mu_0, double *rho_mu_0,  double s_e, double r_e, double gama_pole_h, double rho_pole_h, double gama_equator_h, double rho_equator_h, double *Omega_e, double *Omega_max, int n_of_it, int verbose, bool IsDM);

void updateAB(double s_gp[SDIV + 1], double mu[MDIV + 1], double **Omega_diff, struct evolution_variables *ev, struct stellar_properties *props, struct DiffRotParams *DiffRotParam, double r_ratio, double *omega_mu_0, double *rho_mu_0, double s_e, double r_e, double gama_pole_h, double rho_pole_h, double gama_equator_h, double rho_equator_h, int verbose, bool IsDM);

void compute_velocity_energy_pressure(double *s_gp, double *mu, struct evolution_variables *ev, struct stellar_properties *star_props, struct EOS *eos, struct DiffRotParams *DiffRotParams, double r_ratio, double r_ratio_other, double r_e, double s_e, double **Omega_diff, double **velocity_sq, double **enthalpy, double **pressure, double **energy, double gama_pole_h, double rho_pole_h, double *sin_theta, int verbose);

int spin(double s_gp[SDIV + 1], double mu[MDIV + 1], struct EOS *eosBM, struct stellar_properties *star_props, struct evolution_variables *ev, double accuracy, double cf, double r_ratio_BM, double *r_e_BM_new, double *Omega_BM, double s_e_BM, struct EOS *eosDM, struct stellar_properties *DM_props, double r_ratio_DM, double *r_e_DM_new, double *Omega_DM, double s_e_DM, int verbose, bool *outOfiter, bool counter, double *Omega_e_BM, struct DiffRotParams *DiffRotBM, double *Omega_e_DM, struct DiffRotParams *DiffRotDM, bool zero);

void print(double r_ratio, double e_center, double Mass, double Mass_0, double R_e, double Omega, double Omega_K, double J, enum EOSType type);

#endif