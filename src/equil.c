/*****************************************************************************
 *	equil.c
 *
 *		The code in this file is a set of procedures written by
 *	Nikolaos Stergioulas. These are the procedures used to integrate
 *	the ev equations for a rapidly rotating neutron star.
 *
 * 	The most important procedures are:
 *
 *	make_grid:	Create the MDIV x SDIV grid.
 *			MDIV = iber of divisions of variable mu=cos theta
 *			SDIV = iber of divisions of radial variable s
 *	load_eos:	Load the equation of state file
 *	make_center:	Calculate the central pressure and enthalpy
 *	sphere:		Compute the metric of a spherical star
 *	TOV:		Integrates the Tolman-Oppenheimer-Volkoff
 *			equations for spherically symmetric star
 *	spin:		Integrates the equations for a rapidly rotating
 *			neutron star with oblateness = r_ratio =
 *				radius of pole/radius of equator
 *	mass_radius:	Calculates the gravitational mass and equatorial
 *			radius of the rotating star, along with other
 *			equilibrium quantities.
 *
 ******************************************************************************/

#include "equil.h"

#include <fenv.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "consts.h"
#include "equil_util.h"
#include "nrutil.h"
#include "output.h"

/*******************************************************************/
/* Allocate evolution variables for a grid SDIV, MDIV	      	   */
/* 								                                   */
/*       							                               */
/*******************************************************************/

void allocate_evolution_variables(struct evolution_variables *ev) {
    ev->rho = dmatrix(1, SDIV, 1, MDIV);
    ev->gama = dmatrix(1, SDIV, 1, MDIV);
    ev->alpha = dmatrix(1, SDIV, 1, MDIV);
    ev->omega = dmatrix(1, SDIV, 1, MDIV);
    ev->energyBM = dmatrix(1, SDIV, 1, MDIV);
    ev->pressureBM = dmatrix(1, SDIV, 1, MDIV);
    ev->enthalpyBM = dmatrix(1, SDIV, 1, MDIV);
    ev->velocity_sqBM = dmatrix(1, SDIV, 1, MDIV);
    ev->Omega_diffBM = dmatrix(1, SDIV, 1, MDIV);
    ev->energyDM = dmatrix(1, SDIV, 1, MDIV);
    ev->pressureDM = dmatrix(1, SDIV, 1, MDIV);
    ev->enthalpyDM = dmatrix(1, SDIV, 1, MDIV);
    ev->velocity_sqDM = dmatrix(1, SDIV, 1, MDIV);
    ev->Omega_diffDM = dmatrix(1, SDIV, 1, MDIV);
}

void deallocate_evolution_variables(struct evolution_variables *ev) {
    free_dmatrix(ev->rho, 1, SDIV, 1, MDIV);
    free_dmatrix(ev->gama, 1, SDIV, 1, MDIV);
    free_dmatrix(ev->alpha, 1, SDIV, 1, MDIV);
    free_dmatrix(ev->omega, 1, SDIV, 1, MDIV);
    free_dmatrix(ev->energyBM, 1, SDIV, 1, MDIV);
    free_dmatrix(ev->pressureBM, 1, SDIV, 1, MDIV);
    free_dmatrix(ev->enthalpyBM, 1, SDIV, 1, MDIV);
    free_dmatrix(ev->velocity_sqBM, 1, SDIV, 1, MDIV);
    free_dmatrix(ev->Omega_diffBM, 1, SDIV, 1, MDIV);
    free_dmatrix(ev->energyDM, 1, SDIV, 1, MDIV);
    free_dmatrix(ev->pressureDM, 1, SDIV, 1, MDIV);
    free_dmatrix(ev->enthalpyDM, 1, SDIV, 1, MDIV);
    free_dmatrix(ev->velocity_sqDM, 1, SDIV, 1, MDIV);
    free_dmatrix(ev->Omega_diffDM, 1, SDIV, 1, MDIV);
}
/*******************************************************************/
/* Create computational grid.                                      */
/* Points in the mu-direction are stored in the array mu[i].       */
/* Points in the s-direction are stored in the array s_gp[j].      */
/*******************************************************************/

void make_grid(double s_gp[SDIV + 1], double mu[MDIV + 1]) {
    for (int s = 1; s <= SDIV; s++)
        s_gp[s] = SMAX * (s - 1.0) / (SDIV - 1.0);
    for (int m = 1; m <= MDIV; m++)
        mu[m] = (m - 1.0) / (MDIV - 1.0);
}

/*************************************************************************/
/* EOS management                                                        */
/*************************************************************************/

void load_eos(struct EOS *eos) {
    int i;
    double p, rho, h, n0;
    FILE *f_eos;
    if ((f_eos = fopen(eos->filepath, "r")) == NULL) {
        printf("Cannot open file:  %s\n", eos->filepath);
        exit(1);
    }

    int n = fscanf(f_eos, "%d", &(eos->n_tab));
    if (n != 1) {
        printf("The first line must be one integer!");
        exit(1);
    }

    /* Here it reads the EoS, making the quantities dimensionless */
    for (i = 1; i <= eos->n_tab; i++) {
        n = fscanf(f_eos, "%lf %lf %lf %lf\n", &rho, &p, &h, &n0);
        if (n != 4) {
            printf("Something wrong with the EoS table!");
            printf("In %s\n", eos->filepath);
            printf("%d values out of 4 were found", n);
            exit(1);
        }
        /*To get energy density, one needs to multiply by c squared*/
        eos->log_e_tab[i] = log10(rho * C * C * KSCALE);
        eos->log_p_tab[i] = log10(p * KSCALE);
        eos->log_h_tab[i] = log10(h / (C * C));
        eos->log_n0_tab[i] = log10(n0);
    }
}

/*******************************************************************/
double e_of_rho0(double rho0, double Gamma_P) {
    return (pow(rho0, Gamma_P) / (Gamma_P - 1.0) + rho0);
}

/*C*/
/*******************************************************************/
double e_at_p(double pp, struct EOS *eos, int *n_nearest_pt) {
    if (eos->type == TAB)
        return pow(10.0, interp(eos->log_p_tab, eos->log_e_tab, eos->n_tab, log10(pp), n_nearest_pt));
    else
        return pp / (eos->Gamma_P - 1.0) + pow(pp, 1.0 / eos->Gamma_P);
}

/*C*/
/*******************************************************************/

double p_at_e(double ee, struct EOS *eos, int *n_nearest_pt) {
    return pow(10.0, interp(eos->log_e_tab, eos->log_p_tab, eos->n_tab, log10(ee), n_nearest_pt));
}

/*C*/
/*******************************************************************/

double p_at_h(double hh, struct EOS *eos, int *n_nearest_pt) {
    return pow(10.0, interp(eos->log_h_tab, eos->log_p_tab, eos->n_tab, log10(hh), n_nearest_pt));
}

/*C*/
/*******************************************************************/

double h_at_p(double pp, struct EOS *eos, int *n_nearest_pt) {
    return pow(10.0, interp(eos->log_p_tab, eos->log_h_tab, eos->n_tab, log10(pp), n_nearest_pt));
}

/*C*/
/*******************************************************************/

double n0_at_e(double ee, struct EOS *eos, int *n_nearest_pt) {
    return pow(10.0, interp(eos->log_e_tab, eos->log_n0_tab, eos->n_tab, log10(ee), n_nearest_pt));
}

/*C*/
/***************************************************************/
/*******************************************************************/
double e_of_rho0_wrapper(double rho0, double *args) {
    double Gamma_P = args[0];
    double e_target = args[1];
    return e_of_rho0(rho0, Gamma_P) - e_target;
}
void make_center(struct EOS *eos, double e_center, double *p_center, double *h_center) {
    int n_nearest;
    double rho0_center;

    n_nearest = eos->n_tab / 2;

    if (eos->type == TAB) {
        (*p_center) = p_at_e(e_center, eos, &n_nearest);
        (*h_center) = h_at_p(*p_center, eos, &n_nearest);
    } else {
        rho0_center = brent_root_finder(e_of_rho0_wrapper, 0.0, e_center, DBL_EPSILON, "rho0_center", 2, eos->Gamma_P, e_center);
        (*p_center) = pow(rho0_center, eos->Gamma_P);
        (*h_center) = log((e_center + (*p_center)) / rho0_center);
    }
}

/*C*/
/***********************************************************************/
/* Computes the gravitational mass, equatorial radius, angular momentum
 *	of the star
 * 	and the velocity of co- and counter-rotating particles
 *	with respect to a ZAMO                                         */
/***********************************************************************/

void mass_radius(double s_gp[SDIV + 1], double mu[MDIV + 1], struct EOS *eosBM, struct evolution_variables *ev, double r_ratio_BM, double r_e_BM, double Omega_BM, double *Mass_BM, double *Mass_0_BM, double *ang_mom_BM, double *R_e_BM, double *v_plus, double *v_minus, double *Omega_K_BM, struct stellar_properties *star_props, struct EOS *eosDM, double r_ratio_DM, double r_e_DM, double Omega_DM, double *Mass_DM, double *Mass_0_DM, double *ang_mom_DM, double *R_e_DM, struct stellar_properties *DM_props, double MDM, double *R_p_BM, double *R_p_DM, double *Omega_K_DM, double *T_BM, double *T_DM, double *W_BM, double *W_DM, bool counter, double *GRV2, double *GRV3) {
    int s, m, n_nearest;
    double **rho_0_BM, **rho_0_DM, **velocity, gama_equator_BM, rho_equator_BM, omega_equator_BM, gama_equator_DM, rho_equator_DM, omega_equator_DM, s1, s_1, d_gama_s, d_rho_s, d_omega_s, sqrt_v, D_m_BM[SDIV + 1] = {0.0}, D_m_DM[SDIV + 1] = {0.0}, D_m_0_BM[SDIV + 1] = {0.0}, D_m_0_DM[SDIV + 1] = {0.0}, D_J_BM[SDIV + 1] = {0.0}, D_J_DM[SDIV + 1] = {0.0}, d_o_e[SDIV + 1], d_g_e[SDIV + 1], d_r_e[SDIV + 1], d_v_e[SDIV + 1], doe, dge, dre, dve, vek, gama_mu_0[SDIV + 1], rho_mu_0[SDIV + 1], omega_mu_0[SDIV + 1], J_BM, J_DM, r_p, s_p, gama_mu_1[SDIV + 1], rho_mu_1[SDIV + 1], D_m_p_BM[SDIV + 1] = {0.0}, D_m_p_DM[SDIV + 1] = {0.0}, D_T_BM[SDIV + 1] = {0.0}, D_T_DM[SDIV + 1] = {0.0};

    double Mass_p_BM, Mass_p_DM, TBM_diff, TDM_diff;
    double r_e_BM_old, r_p_BM, s_p_BM, s_e_BM, r_e_DM_old, r_p_DM, s_p_DM, s_e_DM;

    double r_e = fmax(r_e_BM, r_e_DM);

    // BM
    r_e_BM_old = r_e_BM;
    r_p_BM = r_ratio_BM * r_e_BM;
    s_p_BM = r_p_BM / (r_p_BM + r_e);
    s_e_BM = r_e_BM / (r_e_BM + r_e);

    // DM
    r_e_DM_old = r_e_DM;
    r_p_DM = r_ratio_DM * r_e_DM;
    s_p_DM = r_p_DM / (r_p_DM + r_e);
    s_e_DM = r_e_DM / (r_e_DM + r_e);

    rho_0_BM = dmatrix(1, SDIV, 1, MDIV);
    rho_0_DM = dmatrix(1, SDIV, 1, MDIV);
    velocity = dmatrix(1, SDIV, 1, MDIV);

    for (s = 1; s <= SDIV; s++) {
        gama_mu_0[s] = ev->gama[s][1];
        rho_mu_0[s] = ev->rho[s][1];
        gama_mu_1[s] = ev->gama[s][MDIV];
        rho_mu_1[s] = ev->rho[s][MDIV];
    }

    /* Circumferential radius */
    // BM
    n_nearest = SDIV / 2;
    int n_nearestBM = eosBM->n_tab / 2;
    gama_equator_BM = interp(s_gp, gama_mu_0, SDIV, s_e_BM, &n_nearest);
    rho_equator_BM = interp(s_gp, rho_mu_0, SDIV, s_e_BM, &n_nearest);

    (*R_e_BM) = r_e_BM * exp((gama_equator_BM - rho_equator_BM) / 2.0);
    if (eosBM->type == TAB)
        (*R_e_BM) *= sqrt(KAPPA);

    double gama_pole_BM = interp(s_gp, gama_mu_1, SDIV, s_p_BM, &n_nearest);
    double rho_pole_BM = interp(s_gp, rho_mu_1, SDIV, s_p_BM, &n_nearest);
    *R_p_BM = r_p_BM * exp((gama_pole_BM - rho_pole_BM) / 2.0);
    if (eosBM->type == TAB)
        *R_p_BM *= sqrt(KAPPA);

    // DM
    n_nearest = SDIV / 2;
    int n_nearestDM = eosDM->n_tab / 2;
    gama_equator_DM = interp(s_gp, gama_mu_0, SDIV, s_e_DM, &n_nearest);
    rho_equator_DM = interp(s_gp, rho_mu_0, SDIV, s_e_DM, &n_nearest);

    (*R_e_DM) = r_e_DM * exp((gama_equator_DM - rho_equator_DM) / 2.0);
    if (eosDM->type == TAB)
        (*R_e_DM) *= sqrt(KAPPA);

    double gama_pole_DM = interp(s_gp, gama_mu_1, SDIV, s_p_DM, &n_nearest);
    double rho_pole_DM = interp(s_gp, rho_mu_1, SDIV, s_p_DM, &n_nearest);
    *R_p_DM = r_p_DM * exp((gama_pole_DM - rho_pole_DM) / 2.0);
    if (eosDM->type == TAB)
        *R_p_DM *= sqrt(KAPPA);

    /* Masses and angular momentum */
    (*Mass_BM) = 0.0;
    (*Mass_0_BM) = 0.0;
    (*Mass_DM) = 0.0;
    (*Mass_0_DM) = 0.0;
    Mass_p_BM = 0;
    Mass_p_DM = 0;
    J_BM = 0.0;
    J_DM = 0.0;
    TBM_diff = 0.0;
    TDM_diff = 0.0;

    /* CALCULATE THE REST MASS DENSITY */
    if (eosBM->type == TAB) {
        n_nearest = eosBM->n_tab / 2;
        for (s = 1; s <= SDIV; s++)
            for (m = 1; m <= MDIV; m++) {
                if (ev->energyBM[s][m] > star_props->e_surface)
                    rho_0_BM[s][m] = n0_at_e(ev->energyBM[s][m], eosBM, &n_nearestBM) * MB * KSCALE * SQ(C);
                else
                    rho_0_BM[s][m] = 0.0;
            }
    } else {
        for (s = 1; s <= SDIV; s++)
            for (m = 1; m <= MDIV; m++)
                rho_0_BM[s][m] = (ev->energyBM[s][m] + ev->pressureBM[s][m]) * exp(-ev->enthalpyBM[s][m]);
    }

    if (eosDM->type == TAB) {
        n_nearest = eosDM->n_tab / 2;
        for (s = 1; s <= SDIV; s++)
            for (m = 1; m <= MDIV; m++) {
                if (ev->energyDM[s][m] > DM_props->e_surface)
                    rho_0_DM[s][m] = n0_at_e(ev->energyDM[s][m], eosDM, &n_nearestDM) * MDM * KSCALE * SQ(C);
                else
                    rho_0_DM[s][m] = 0.0;
            }
    } else {
        for (s = 1; s <= SDIV; s++)
            for (m = 1; m <= MDIV; m++)
                rho_0_DM[s][m] = (ev->energyDM[s][m] + ev->pressureDM[s][m]) * exp(-ev->enthalpyDM[s][m]);
    }

    // Simpson's rule weights and overall factor:
    const double weights[3] = {1.0, 4.0, 1.0};
    const double pref_m = 1.0 / (3.0 * (MDIV - 1));
    const double pref_s = SMAX / (3.0 * (SDIV - 1));

    for (int s = 1; s <= SDIV; s++) {
        for (int m = 1; m <= MDIV - 2; m += 2) {
            double acc_DmBM = 0, acc_DmDM = 0;
            double acc_Dm0BM = 0, acc_Dm0DM = 0, acc_DmpBM = 0, acc_DmpDM = 0;
            double acc_JBM = 0, acc_JDM = 0, acc_TBM = 0, acc_TDM = 0;

            for (int i = 0; i < 3; i++) {
                int mi = m + i;

                // Precompute common expressions for indices m, m+1, and m+2:
                double exp_alpha_gama_m = exp(2.0 * ev->alpha[s][mi] + ev->gama[s][mi]);
                double exp_alpha_gama_rho_m = exp(2.0 * ev->alpha[s][mi] + (ev->gama[s][mi] - ev->rho[s][mi]) / 2.0);
                double exp_alpha_gama_m_rho_m = exp(2.0 * ev->alpha[s][mi] + ev->gama[s][mi] - ev->rho[s][mi]);
                double sqrt_velocity_sqBM_m = sqrt(ev->velocity_sqBM[s][mi]);
                double sqrt_velocity_sqDM_m = sqrt(ev->velocity_sqDM[s][mi]);
                double sqrt_mu_m = sqrt(1.0 - mu[mi] * mu[mi]);

                // --- Compute contributions for gravitational mass:
                acc_DmBM += weights[i] * exp_alpha_gama_m * (((ev->energyBM[s][mi] + ev->pressureBM[s][mi]) / (1.0 - ev->velocity_sqBM[s][mi])) * (1.0 + ev->velocity_sqBM[s][mi] + (2.0 * s_gp[s] * sqrt_velocity_sqBM_m / (1.0 - s_gp[s])) * sqrt_mu_m * r_e * ev->omega[s][mi] * exp(-ev->rho[s][mi])) + 2.0 * ev->pressureBM[s][mi]);

                acc_DmDM += weights[i] * exp_alpha_gama_m * (((ev->energyDM[s][mi] + ev->pressureDM[s][mi]) / (1.0 - ev->velocity_sqDM[s][mi])) * (1.0 + ev->velocity_sqDM[s][mi] + (2.0 * s_gp[s] * sqrt_velocity_sqDM_m / (1.0 - s_gp[s])) * sqrt(1.0 - mu[mi] * mu[mi]) * r_e * ev->omega[s][mi] * exp(-ev->rho[s][mi])) + 2.0 * ev->pressureDM[s][mi]);

                // --- Compute contributions for the rest mass:
                acc_Dm0BM += weights[i] * exp_alpha_gama_rho_m * (rho_0_BM[s][mi] / sqrt(1.0 - ev->velocity_sqBM[s][mi]));
                acc_Dm0DM += weights[i] * exp_alpha_gama_rho_m * (rho_0_DM[s][mi] / sqrt(1.0 - ev->velocity_sqDM[s][mi]));

                // --- Compute contributions for the proper mass:
                acc_DmpBM += weights[i] * exp_alpha_gama_rho_m * (ev->energyBM[s][mi] / sqrt(1.0 - ev->velocity_sqBM[s][mi]));
                acc_DmpDM += weights[i] * exp_alpha_gama_rho_m * (ev->energyDM[s][mi] / sqrt(1.0 - ev->velocity_sqDM[s][mi]));

                // --- Compute contributions for angular momentum:
                acc_JBM += weights[i] * sqrt_mu_m * exp_alpha_gama_m_rho_m * (ev->energyBM[s][mi] + ev->pressureBM[s][mi]) * sqrt_velocity_sqBM_m / (1.0 - ev->velocity_sqBM[s][mi]);

                acc_JDM += weights[i] * sqrt_mu_m * exp_alpha_gama_m_rho_m * (ev->energyDM[s][mi] + ev->pressureDM[s][mi]) * sqrt_velocity_sqDM_m / (1.0 - ev->velocity_sqDM[s][mi]);

                // --- Compute rotational kinetic energy:
                acc_TBM += weights[i] * sqrt_mu_m * exp_alpha_gama_m_rho_m * (ev->energyBM[s][mi] + ev->pressureBM[s][mi]) * sqrt_velocity_sqBM_m * ev->Omega_diffBM[s][mi] / (1 - ev->velocity_sqBM[s][mi]);

                acc_TDM += weights[i] * sqrt_mu_m * exp_alpha_gama_m_rho_m * (ev->energyDM[s][mi] + ev->pressureDM[s][mi]) * sqrt_velocity_sqDM_m * ev->Omega_diffDM[s][mi] / (1 - ev->velocity_sqDM[s][mi]);
            }

            D_m_BM[s] += acc_DmBM;
            D_m_DM[s] += acc_DmDM;
            D_m_0_BM[s] += acc_Dm0BM;
            D_m_0_DM[s] += acc_Dm0DM;
            D_m_p_BM[s] += acc_DmpBM;
            D_m_p_DM[s] += acc_DmpDM;
            D_J_BM[s] += acc_JBM;
            D_J_DM[s] += acc_JDM;
            D_T_BM[s] += acc_TBM;
            D_T_DM[s] += acc_TDM;
        }
    }

    for (s = 1; s <= SDIV - 2; s += 2) {
        for (int i = 0; i < 3; i++) {
            int si = s + i;

            double inv = 1.0 / (1.0 - s_gp[si]);
            double tmp1 = SQ(s_gp[s]) * SQ(SQ(inv));                    // s^2 / (1 - s)^4
            double tmp4 = (SQ(s_gp[s]) * s_gp[s]) * SQ(SQ(inv)) * inv;  // s^3 / (1 - s)^5

            *Mass_BM += weights[i] * tmp1 * D_m_BM[si];
            *Mass_DM += weights[i] * tmp1 * D_m_DM[si];
            *Mass_0_BM += weights[i] * tmp1 * D_m_0_BM[si];
            *Mass_0_DM += weights[i] * tmp1 * D_m_0_DM[si];
            Mass_p_BM += weights[i] * tmp1 * D_m_p_BM[si];
            Mass_p_DM += weights[i] * tmp1 * D_m_p_DM[si];
            J_BM += weights[i] * tmp4 * D_J_BM[si];
            J_DM += weights[i] * tmp4 * D_J_DM[si];
            TBM_diff += weights[i] * tmp4 * D_T_BM[si];
            TDM_diff += weights[i] * tmp4 * D_T_DM[si];
        }
    }

    // Apply integration prefactors
    *Mass_BM *= pref_s * pref_m;
    *Mass_DM *= pref_s * pref_m;
    *Mass_0_BM *= pref_s * pref_m;
    *Mass_0_DM *= pref_s * pref_m;
    Mass_p_BM *= pref_s * pref_m;
    Mass_p_DM *= pref_s * pref_m;
    J_BM *= pref_s * pref_m;
    J_DM *= pref_s * pref_m;
    TBM_diff *= pref_s * pref_m;
    TDM_diff *= pref_s * pref_m;

    double re3 = r_e * r_e * r_e;
    (*Mass_BM) *= 4 * PI * re3;
    (*Mass_0_BM) *= 4 * PI * re3;
    Mass_p_BM *= 4 * PI * re3;
    if (eosBM->type == TAB) {
        (*Mass_BM) *= sqrt(KAPPA) * C * C / G;
        (*Mass_0_BM) *= sqrt(KAPPA) * C * C / G;
        Mass_p_BM *= sqrt(KAPPA) * C * C / G;
    }

    (*Mass_DM) *= 4 * PI * re3;
    (*Mass_0_DM) *= 4 * PI * re3;
    Mass_p_DM *= 4 * PI * re3;
    if (eosDM->type == TAB) {
        (*Mass_DM) *= sqrt(KAPPA) * C * C / G;
        (*Mass_0_DM) *= sqrt(KAPPA) * C * C / G;
        Mass_p_DM *= sqrt(KAPPA) * C * C / G;
    }

    double re4 = re3 * r_e;
    if (r_ratio_BM == 1.0 && r_ratio_DM == 1.0) {
        J_BM = 0.0;
    } else {
        J_BM *= 4 * PI * re4;
        TBM_diff *= 2.0 * PI * re4;
        if (eosBM->type == TAB) {
            J_BM *= KAPPA * C * C * C / G;
            TBM_diff *= sqrt(KAPPA) * C * C * C * C / G;
        }
    }

    if (r_ratio_DM == 1.0 && r_ratio_BM == 1.0) {
        J_DM = 0.0;
    } else {
        J_DM *= 4 * PI * re4;
        TDM_diff *= 2.0 * PI * re4;
        if (eosDM->type == TAB) {
            J_DM *= KAPPA * C * C * C / G;
            TDM_diff *= sqrt(KAPPA) * C * C * C * C / G;
        }
    }

    (*ang_mom_BM) = J_BM;
    (*ang_mom_DM) = J_DM;

    if (star_props->RotType == UNIFORM) {
        *T_BM = 0.5 * J_BM * Omega_BM;
    } else {
        *T_BM = TBM_diff;
    }

    if (DM_props->RotType == UNIFORM) {
        *T_DM = 0.5 * J_DM * Omega_DM;
    } else {
        *T_DM = TDM_diff;
    }

    *W_BM = Mass_p_BM * C * C - *Mass_BM * C * C + *T_BM;
    *W_DM = Mass_p_DM * C * C - *Mass_DM * C * C + *T_DM;

    /* Compute the velocities of co-rotating and counter-rotating particles with respect to a ZAMO 	*/

    for (s = 1 + (SDIV - 1) / 2; s <= SDIV; s++) {
        s1 = s_gp[s] * (1.0 - s_gp[s]);
        s_1 = 1.0 - s_gp[s];

        d_gama_s = deriv_s(ev->gama, s, 1);
        d_rho_s = deriv_s(ev->rho, s, 1);
        d_omega_s = deriv_s(ev->omega, s, 1);

        sqrt_v = exp(-2.0 * ev->rho[s][1]) * r_e * r_e * pow(s_gp[s], 4.0) * pow(d_omega_s, 2.0) + 2 * s1 * (d_gama_s + d_rho_s) + s1 * s1 * (d_gama_s * d_gama_s - d_rho_s * d_rho_s);

        if (sqrt_v > 0.0) {
            sqrt_v = sqrt(sqrt_v);
        } else {
            sqrt_v = 0.0;
        }

        v_plus[s] = (exp(-ev->rho[s][1]) * r_e * s_gp[s] * s_gp[s] * d_omega_s + sqrt_v) /
                    (2.0 + s1 * (d_gama_s - d_rho_s));

        v_minus[s] = (exp(-ev->rho[s][1]) * r_e * s_gp[s] * s_gp[s] * d_omega_s - sqrt_v) /
                     (2.0 + s1 * (d_gama_s - d_rho_s));
    }

    /* Kepler angular velocity */

    for (s = 1; s <= SDIV; s++) {
        d_o_e[s] = deriv_s(ev->omega, s, 1);
        d_g_e[s] = deriv_s(ev->gama, s, 1);
        d_r_e[s] = deriv_s(ev->rho, s, 1);
        d_v_e[s] = deriv_s(velocity, s, 1);
        /* Value of omega on the equatorial plane*/
        omega_mu_0[s] = ev->omega[s][1];
    }

    n_nearest = SDIV / 2;

    // General formula:
    // Eq. 23b of https://ui.adsabs.harvard.edu/abs/1986ApJ...304..115F/abstract
    double sv = s_e_BM;
    doe = interp(s_gp, d_o_e, SDIV, sv, &n_nearest);
    dge = interp(s_gp, d_g_e, SDIV, sv, &n_nearest);
    dre = interp(s_gp, d_r_e, SDIV, sv, &n_nearest);
    dve = interp(s_gp, d_v_e, SDIV, sv, &n_nearest);
    double rho_equator = interp(s_gp, rho_mu_0, SDIV, sv, &n_nearest);
    double omega_equator = interp(s_gp, omega_mu_0, SDIV, sv, &n_nearest);
    double tmp1 = (sv * sv * doe / (2.0 + (dge - dre) * sv * (1 - sv))) * r_e * exp(-rho_equator);
    double tmp2 = sv * (1 - sv) * (dge + dre) / (2.0 + (dge - dre) * sv * (1 - sv));
    vek = tmp1 + sqrt(tmp2 + SQ(tmp1));

    // Eq. 23a of https://ui.adsabs.harvard.edu/abs/1986ApJ...304..115F/abstract
    (*Omega_K_BM) = (C / sqrt(KAPPA)) * (omega_equator + vek * exp(rho_equator) / r_e_BM);

    sv = s_e_DM;
    doe = interp(s_gp, d_o_e, SDIV, sv, &n_nearest);
    dge = interp(s_gp, d_g_e, SDIV, sv, &n_nearest);
    dre = interp(s_gp, d_r_e, SDIV, sv, &n_nearest);
    dve = interp(s_gp, d_v_e, SDIV, sv, &n_nearest);
    rho_equator = interp(s_gp, rho_mu_0, SDIV, sv, &n_nearest);
    omega_equator = interp(s_gp, omega_mu_0, SDIV, sv, &n_nearest);
    tmp1 = (sv * sv * doe / (2.0 + (dge - dre) * sv * (1 - sv))) * r_e * exp(-rho_equator);
    tmp2 = sv * (1 - sv) * (dge + dre) / (2.0 + (dge - dre) * sv * (1 - sv));

    // Eq. 23a of https://ui.adsabs.harvard.edu/abs/1986ApJ...304..115F/abstract
    if (counter) {
        vek = tmp1 - sqrt(tmp2 + SQ(tmp1));
    } else {
        vek = tmp1 + sqrt(tmp2 + SQ(tmp1));
    }

    if (r_e_DM > 0) {
        (*Omega_K_DM) = (C / sqrt(KAPPA)) * (omega_equator + vek * exp(rho_equator) / r_e_DM);
    } else {
        *Omega_K_DM = 0.0;
    }
    free_dmatrix(velocity, 1, SDIV, 1, MDIV);
    free_dmatrix(rho_0_BM, 1, SDIV, 1, MDIV);
    free_dmatrix(rho_0_DM, 1, SDIV, 1, MDIV);

    /* GRV2 AND GRV3 */
    compute_virial(ev, s_gp, mu, r_e, r_ratio_BM, GRV2, GRV3);
}

#define SAFE_DERIV_S(field, s, m) ((s) == 1 ? 0.0 : deriv_s(field, s, m))
#define SAFE_DERIV_SS(field, s, m) ((s == 1 || s == 2) ? 0.0 : deriv_ss(field, s, m))
#define SAFE_DERIV_M(field, s, m) ((s) == 1 ? 0.0 : deriv_m(field, s, m))
#define SAFE_DERIV_MM(field, s, m) ((s == 1 || s == 2) ? 0.0 : deriv_s(field, s, m))

void compute_virial(struct evolution_variables *ev, double *s_gp, double *mu, double r_e, double r_ratio, double *grv2, double *grv3) {
    // Accumulators for Simpson integration
    double grv2_virial1 = 0.0, grv2_virial2 = 0.0, grv2_virial3 = 0.0;
    double grv3_virial1 = 0.0, grv3_virial2 = 0.0;

    // Temporary storage for s
    double Dv1_grv2_mu[SDIV + 1] = {0};
    double Dv2_grv2_mu[SDIV + 1] = {0};
    double Dv3_grv2_mu[SDIV + 1] = {0};
    double Dv1_grv3_mu[SDIV + 1] = {0};
    double Dv2_grv3_mu[SDIV + 1] = {0};

    // Build nu and beta metric fields
    double **nu = dmatrix(1, SDIV, 1, MDIV);
    double **beta = dmatrix(1, SDIV, 1, MDIV);
    for (int s = 1; s <= SDIV; ++s) {
        for (int m = 1; m <= MDIV; ++m) {
            nu[s][m] = 2.0 * (ev->gama[s][m] + ev->rho[s][m]);
            beta[s][m] = 2.0 * (ev->gama[s][m] - ev->rho[s][m]);
        }
    }

    // Simpson weights and prefactors
    const double weights[3] = {1.0, 4.0, 1.0};
    const double pref_m = 1.0 / (3.0 * (MDIV - 1));
    const double pref_s = SMAX / (3.0 * (SDIV - 1));

    // First stage: compute Dv arrays for grv2 and grv3
    for (int s = 1; s <= SDIV; s++) {
        double sgp = s_gp[s];
        double oms = 1.0 - sgp;

        for (int m = 1; m <= MDIV - 2; m += 2) {
            double acc1_2 = 0.0, acc2_2 = 0.0, acc3_2 = 0.0;
            double acc1_3 = 0.0, acc2_3 = 0.0;

            for (int i = 0; i < 3; i++) {
                int mi = m + i;
                double omm2 = 1.0 - SQ(mu[mi]) + ((m == MDIV - 2) && (i == 2) ? DBL_EPSILON : 0);
                double sqrt_o2 = sqrt(omm2);

                // All quantities
                double pBM = ev->pressureBM[s][mi];
                double eBM = ev->energyBM[s][mi];
                double v2BM = ev->velocity_sqBM[s][mi];
                double pDM = ev->pressureDM[s][mi];
                double eDM = ev->energyDM[s][mi];
                double v2DM = ev->velocity_sqDM[s][mi];

                double omega = ev->omega[s][mi];
                double OmeBM = ev->Omega_diffBM[s][mi];
                double OmeDM = ev->Omega_diffDM[s][mi];
                double alp_v = ev->alpha[s][mi];
                double nu_v = nu[s][mi];
                double beta_v = beta[s][mi];

                // Derivatives common to both grv2 & grv3
                double d_nu_s = SAFE_DERIV_S(nu, s, mi);
                double d_nu_m = SAFE_DERIV_M(nu, s, mi);
                double d_nu_ss = SAFE_DERIV_SS(nu, s, mi);
                double d_nu_mm = SAFE_DERIV_MM(nu, s, mi);

                double d_b_s = SAFE_DERIV_S(beta, s, mi);
                double d_b_m = SAFE_DERIV_M(beta, s, mi);
                double d_b_ss = SAFE_DERIV_SS(beta, s, mi);
                double d_b_mm = SAFE_DERIV_MM(beta, s, mi);

                double d_o_s = SAFE_DERIV_S(ev->omega, s, mi);
                double d_o_m = SAFE_DERIV_M(ev->omega, s, mi);
                double d_a_s = SAFE_DERIV_S(ev->alpha, s, mi);

                // --- grv2: Dv1 term ---
                double hBM = (eBM + pBM) * v2BM / (1.0 - v2BM);
                double hDM = (eDM + pDM) * v2DM / (1.0 - v2DM);
                acc1_2 += weights[i] * (pBM + hBM + pDM + hDM) * exp(2 * alp_v) / sqrt_o2;

                // --- grv2: Dv2 term ---
                double grad_nu = SQ(oms * d_nu_s) + omm2 * SQ(d_nu_m / sgp);
                double lap_nu = SQ(oms) * d_nu_ss - 2.0 * oms * d_nu_s + oms * d_nu_s / sgp - omm2 * d_nu_mm / SQ(sgp) + mu[mi] * d_nu_m / SQ(sgp);
                acc2_2 += weights[i] * (grad_nu - nu_v * lap_nu) / sqrt_o2;

                // --- grv2: Dv3 term ---
                double A_BM = (OmeBM - omega) * (eBM + pBM) / (1.0 - v2BM);
                double A_DM = (OmeDM - omega) * (eDM + pDM) / (1.0 - v2DM);
                double grad_o = SQ(oms) * oms / sgp * d_o_s - mu[mi] * SQ(oms) * d_o_m / SQ(sgp) + SQ(SQ(oms)) * d_o_s * (d_b_s + d_nu_s) + SQ(oms / sgp) * omm2 * d_o_m * (d_b_m + d_nu_m);
                acc3_2 += weights[i] * (16.0 * PI * SQ(r_e) * exp(2 * alp_v) * (A_BM + A_DM) + grad_o) * sqrt_o2 * exp(2.0 * (beta_v - nu_v)) * omega;

                // --- grv3: Dv1 term (3 * pressure + enthalpy factor) ---
                double t1 = 3.0 * pBM + (eBM + pBM) * v2BM / (1.0 - v2BM) + 3.0 * pDM + (eDM + pDM) * v2DM / (1.0 - v2DM);
                acc1_3 += weights[i] * t1 * exp(2 * alp_v + beta_v);

                // --- grv3: Dv2 term (complex metric derivatives) ---
                double term1 = SQ(SQ(oms)) * d_nu_ss / SQ(r_e) - SQ(oms) * omm2 * d_nu_mm / SQ(sgp * r_e) - 2.0 * SQ(oms) * oms * d_nu_s / SQ(r_e) + 2.0 * SQ(oms) * oms * d_nu_s / sgp / SQ(r_e);
                double term2 = SQ(SQ(oms)) * d_nu_s * d_b_s / SQ(r_e) + SQ(oms) * omm2 * d_nu_m * d_b_m / SQ(sgp * r_e);
                double term3 = SQ(SQ(oms)) * d_b_ss / SQ(r_e) - SQ(oms) * omm2 * d_b_mm / SQ(sgp * r_e) - 2.0 * SQ(oms) * oms * d_b_s / SQ(r_e) + 2.0 * SQ(oms) * oms * d_b_s / sgp / SQ(r_e);
                double term4 = SQ(SQ(oms)) * d_b_s * d_b_s / SQ(r_e) + SQ(oms) * omm2 * d_b_m * d_b_m / SQ(sgp * r_e);
                acc2_3 += weights[i] * exp(beta_v) * (-nu_v * (term1 + term2) + 0.5 * alp_v * (term3 + term4) - SQ(oms) * oms * mu[mi] * d_a_s * (1 - exp(2 * alp_v - 2 * beta_v)) / (2.0 * SQ(r_e) * sgp * sqrt_o2) + SQ(oms) * oms * mu[mi] * d_b_s * (1 - exp(2 * alp_v - 2 * beta_v)) / (4.0 * SQ(r_e) * sgp * sqrt_o2) - 6.0 * PI * omega * exp(4.0 * (beta_v - nu_v)) * ((eBM + pBM) * (OmeBM - omega) / (1.0 - v2BM) + (eDM + pDM) * (OmeDM - omega) / (1.0 - v2DM)) * SQ(r_e * sgp / oms) * omm2);
            }

            Dv1_grv2_mu[s] += acc1_2;
            Dv2_grv2_mu[s] += acc2_2;
            Dv3_grv2_mu[s] += acc3_2;
            Dv1_grv3_mu[s] += acc1_3;
            Dv2_grv3_mu[s] += acc2_3;
        }
    }

    // Second stage: Simpson over s-index with inner i-loop
    for (int s = 1; s <= SDIV - 2; s += 2) {
        for (int i = 0; i < 3; i++) {
            int si = s + i;
            double inv = 1.0 / (1.0 - s_gp[si]);
            double s0 = s_gp[si];
            double w1 = s0 * SQ(inv) * inv;                      // inv^3
            double w2 = s0 * inv;                                // inv^1
            double w3 = s0 * s0 * s0 * SQ(inv) * SQ(inv) * inv;  // inv^5
            double wg = s0 * s0 * SQ(inv) * SQ(inv);             // inv^4

            grv2_virial1 += weights[i] * w1 * Dv1_grv2_mu[si];
            grv2_virial2 += weights[i] * w2 * Dv2_grv2_mu[si];
            grv2_virial3 += weights[i] * w3 * Dv3_grv2_mu[si];
            grv3_virial1 += weights[i] * wg * Dv1_grv3_mu[si];
            grv3_virial2 += weights[i] * wg * Dv2_grv3_mu[si];
        }
    }

    // Apply prefactors
    grv2_virial1 *= 2.0 * SQ(r_e) * pref_s * pref_m;
    grv2_virial2 *= pref_s * pref_m;
    grv2_virial3 *= -1.5 * SQ(r_e) * pref_s * pref_m;

    grv3_virial1 *= 2.0 * SQ(r_e) * r_e * pref_s * pref_m;
    grv3_virial2 *= 2.0 * SQ(r_e) * r_e * pref_s * pref_m;

    // printf("grv2 : virial_1 = %g virial_2= %g virial_3 = %g\n", grv2_virial1, grv2_virial2, grv2_virial3);
    // printf("grv3 : virial_1 = %g virial_2= %g\n", grv3_virial1, grv3_virial2);

    // Cleanup
    free_dmatrix(nu, 1, SDIV, 1, MDIV);
    free_dmatrix(beta, 1, SDIV, 1, MDIV);

    *grv2 = fabs(8.0 * PI * grv2_virial1 / (grv2_virial2 + grv2_virial3));
    *grv3 = fabs(4.0 * PI * grv3_virial1 / grv3_virial2);
}

/*C*/
/**************************************************************************/

double dm_dr_is(double r_is, double r, double m, double p, struct stellar_properties *star_props, struct EOS *eos, int *n_nearest_pt) {
    double dmdr,
        e_d;

    if (p < star_props->p_surface)
        e_d = 0.0;
    else
        e_d = e_at_p(p, eos, n_nearest_pt);

    if (r_is < RMIN)
        dmdr = 4.0 * PI * star_props->e_center * r * r * (1.0 + 4.0 * PI * star_props->e_center * r * r / 3.0);
    else
        dmdr = 4.0 * PI * e_d * r * r * r * sqrt(1.0 - 2.0 * m / r) / r_is;

    return dmdr;
}

/*C*/
/**************************************************************************/
double dp_dr_is(double r_is, double r, double m, double p, double ptot, struct stellar_properties *star_props, struct stellar_properties *DM_props, struct EOS *eos, int *n_nearest_pt) {
    double dpdr, e_d;

    if (p < star_props->p_surface)
        e_d = 0.0;
    else
        e_d = e_at_p(p, eos, n_nearest_pt);

    if (r_is < RMIN) {
        dpdr = -4.0 * PI * (star_props->e_center + p) * (star_props->e_center + DM_props->e_center + 3.0 * ptot) * r * (1.0 + 4.0 * (star_props->e_center + DM_props->e_center) * r * r / 3.0) / 3.0;
    } else
        dpdr = -(e_d + p) * (m + 4.0 * PI * r * r * r * ptot) / (r * r_is * sqrt(1.0 - 2.0 * m / r));
    return dpdr;
}

/**************************************************************************/
double dr_dr_is(double r_is, double r, double m) {
    double drdris;

    if (r_is < RMIN)
        drdris = 1.0;
    else
        drdris = (r / r_is) * sqrt(1.0 - 2.0 * m / r);

    return drdris;
}

/*C*/
/************************************************************************/
void rk4_step(double hdiv, double r_is, double r, double m, double p, double m_DM, double p_DM, bool skip_BM, bool skip_DM, struct stellar_properties *star_props, struct stellar_properties *DM_props, struct EOS *eosBM, struct EOS *eosDM, struct rk_coeff *rkc_old, struct rk_coeff *rkc_new) {
    // Set to zero in case we skip one fluid
    rkc_new->b = 0;
    rkc_new->c = 0;
    rkc_new->b_DM = 0;
    rkc_new->c_DM = 0;

    // Euler step with old coefficients
    double r_is_next = r_is + hdiv;
    double r_next = r + hdiv * rkc_old->a;
    double m_next = m + hdiv * rkc_old->b;
    double p_next = p + hdiv * rkc_old->c;
    double m_DM_next = m_DM + hdiv * rkc_old->b_DM;
    double p_DM_next = p_DM + hdiv * rkc_old->c_DM;

    int n_nearest_BM = eosBM->n_tab / 2;
    int n_nearest_DM = eosDM->n_tab / 2;

    double mtot_next = m_next + m_DM_next;
    double ptot_next = p_next + p_DM_next;

    // Compute new coefficients
    rkc_new->a = dr_dr_is(r_is_next, r_next, mtot_next);
    if (!skip_BM) {
        rkc_new->b = dm_dr_is(r_is_next, r_next, mtot_next, p_next, star_props, eosBM, &n_nearest_BM);
        rkc_new->c = dp_dr_is(r_is_next, r_next, mtot_next, p_next, ptot_next, star_props, DM_props, eosBM, &n_nearest_BM);
    }
    if (!skip_DM) {
        rkc_new->b_DM = dm_dr_is(r_is_next, r_next, mtot_next, p_DM_next, DM_props, eosDM, &n_nearest_DM);
        rkc_new->c_DM = dp_dr_is(r_is_next, r_next, mtot_next, p_DM_next, ptot_next, DM_props, star_props, eosDM, &n_nearest_DM);
    }
}

void TOV(int i_check, struct stellar_properties *star_props, struct EOS *eosBM, double r_is_gp[RDIV + 1], double lambda_gp[RDIV + 1], double nu_gp[RDIV + 1], double *r_is_final_big, double *r_final_big, double *m_final, struct stellar_properties *DM_props, struct EOS *eosDM, double *r_is_final_small, double *r_final_small, double *m_DM_final, bool *is_DM) {
    int i = 2, n_nearest_BM, n_nearest_DM;
    double r, r_is, r_is_est, r_is_check, dr_is_save, h, nu_s, hh, k_rescale;
    double e_d, p, m, r_gp[RDIV + 1], m_gp[RDIV + 1], e_d_gp[RDIV + 1], p_d_gp[RDIV + 1];
    double e_d_DM, p_DM, m_DM, r_gp_DM[RDIV + 1], m_gp_DM[RDIV + 1], e_d_gp_DM[RDIV + 1], p_d_gp_DM[RDIV + 1];

    struct rk_coeff rkc[4];

    if (i_check == 1) {
        if ((eosBM->type == TAB) || (eosDM->type == TAB))
            r_is_est = 1.5e6 / sqrt(KAPPA);
        else
            r_is_est = 2.0 * sqrt(eosBM->Gamma_P / (4.0 * PI * (eosBM->Gamma_P - 1.0))) *
                       pow(star_props->e_center, (eosBM->Gamma_P - 2.0) / 2.0);

        h = r_is_est / 10000;
    } else {
        r_is_est = (*r_is_final_big);
        h = r_is_est / 100000;
        dr_is_save = (*r_is_final_big) / RDIV;
        r_is_check = dr_is_save;
    }

    r_is = 0.0; /* initial isotropic radius */
    r = 0.0;    /* initial radius */
    r_is_gp[1] = 0.0;
    r_gp[1] = 0.0;
    lambda_gp[1] = 0.0;

    /* BM */
    m = 0.0;                  /* initial mass */
    p = star_props->p_center; /* initial pressure */

    m_gp[1] = 0.0;
    e_d_gp[1] = star_props->e_center;
    p_d_gp[1] = star_props->p_center;

    /* DM */
    m_DM = 0.0;                /* initial mass */
    p_DM = DM_props->p_center; /* initial pressure */

    m_gp_DM[1] = 0.0;
    e_d_gp_DM[1] = DM_props->e_center;
    p_d_gp_DM[1] = DM_props->p_center;

    n_nearest_BM = eosBM->n_tab / 2;
    n_nearest_DM = eosDM->n_tab / 2;

    bool skip_BM = false,
         skip_DM = false,
         done = false;

    if (DM_props->e_center == 0) {
        skip_DM = true;
    }
    // while (p >= star_props->p_surface || p_DM >= DM_props->p_surface) {
    while (p >= DBL_EPSILON || p_DM >= DBL_EPSILON) {
        if (p < star_props->p_surface) {
            // If below surface pressure, skip integration
            skip_BM = true;
            p = 0;
            e_d = 0;
        } else {
            e_d = e_at_p(p, eosBM, &n_nearest_BM);
        }
        if (p_DM < DM_props->p_surface || skip_DM) {
            // If below surface pressure, skip integration
            skip_DM = true;
            p_DM = 0;
            e_d_DM = 0;
        } else {
            e_d_DM = e_at_p(p_DM, eosDM, &n_nearest_DM);
        }

        if ((i_check == 3) && (r_is > r_is_check) && (i <= RDIV)) {
            r_is_gp[i] = r_is;
            r_gp[i] = r;
            m_gp[i] = m;
            m_gp_DM[i] = m_DM;
            e_d_gp[i] = e_d;
            e_d_gp_DM[i] = e_d_DM;
            p_d_gp[i] = p;
            p_d_gp_DM[i] = p_DM;
            i++;
            r_is_check += dr_is_save;
        }

        // IF (NOT done) AND (skipbm XOR skipdm)
        if (!done && (skip_BM ^ skip_DM)) {
            done = true;
            if (skip_DM) {
                *is_DM = true;
            }
            (*r_is_final_small) = r_is;
            (*r_final_small) = r;
        }
        (*r_is_final_big) = r_is;
        (*r_final_big) = r;
        (*m_final) = m;
        (*m_DM_final) = m_DM;

        rk4_step(0, r_is, r, m, p, m_DM, p_DM, skip_BM, skip_DM, star_props, DM_props, eosBM, eosDM, &rkc[0], &rkc[0]);
        rk4_step(h / 2.0, r_is, r, m, p, m_DM, p_DM, skip_BM, skip_DM, star_props, DM_props, eosBM, eosDM, &rkc[0], &rkc[1]);
        rk4_step(h / 2.0, r_is, r, m, p, m_DM, p_DM, skip_BM, skip_DM, star_props, DM_props, eosBM, eosDM, &rkc[1], &rkc[2]);
        rk4_step(h, r_is, r, m, p, m_DM, p_DM, skip_BM, skip_DM, star_props, DM_props, eosBM, eosDM, &rkc[2], &rkc[3]);

        r += (h / 6.0) * (rkc[0].a + 2 * rkc[1].a + 2 * rkc[2].a + rkc[3].a);
        m += (h / 6.0) * (rkc[0].b + 2 * rkc[1].b + 2 * rkc[2].b + rkc[3].b);
        p += (h / 6.0) * (rkc[0].c + 2 * rkc[1].c + 2 * rkc[2].c + rkc[3].c);
        m_DM += (h / 6.0) * (rkc[0].b_DM + 2 * rkc[1].b_DM + 2 * rkc[2].b_DM + rkc[3].b_DM);
        p_DM += (h / 6.0) * (rkc[0].c_DM + 2 * rkc[1].c_DM + 2 * rkc[2].c_DM + rkc[3].c_DM);

        r_is += h;
    }

    r_is_gp[RDIV] = (*r_is_final_big);
    r_gp[RDIV] = (*r_final_big);
    m_gp[RDIV] = (*m_final);
    m_gp_DM[RDIV] = (*m_DM_final);

    /* Rescale r_is and compute lambda */
    double Ptot[RDIV + 1], Etot[RDIV + 1];
    for (i = 1; i < RDIV + 1; i++) {
        Ptot[i] = p_d_gp[i] + p_d_gp_DM[i];
        Etot[i] = e_d_gp[i] + e_d_gp_DM[i];
    }

    if (i_check == 3) {
        // Scale so that at infinty r = r_is
        k_rescale = 0.5 * ((*r_final_big) / (*r_is_final_big)) * (1.0 - (*m_final + *m_DM_final) / (*r_final_big) + sqrt(1.0 - 2.0 * (*m_final + *m_DM_final) / (*r_final_big)));

        (*r_is_final_big) *= k_rescale;
        (*r_is_final_small) *= k_rescale;

        nu_s = log((1.0 - (*m_final + *m_DM_final) / (2.0 * (*r_is_final_big))) / (1.0 + (*m_final + *m_DM_final) / (2.0 * (*r_is_final_big))));

        for (i = 1; i <= RDIV; i++) {
            r_is_gp[i] *= k_rescale;

            if (i == 1)
                lambda_gp[1] = log(1.0 / k_rescale);
            else
                lambda_gp[i] = log(r_gp[i] / r_is_gp[i]);

            if (Etot[i] < star_props->e_surface + DM_props->e_surface) {
                hh = 0.0;
            } else {
                double sum = 0.0, /* final sum */
                    h0, h1, hph, hdh, hmh;

                /* https://en.wikipedia.org/wiki/Simpson%27s_rule#Composite_Simpson's_rule_for_irregularly_spaced_data */
                for (int j = RDIV; j >= i + 2; j -= 2) {
                    h0 = -(Ptot[j] - Ptot[j - 1]);
                    h1 = -(Ptot[j - 1] - Ptot[j - 2]);

                    hph = h1 + h0;
                    hdh = h1 / h0;
                    hmh = h1 * h0;

                    sum += hph / 6. * ((2. - hdh) / (Ptot[j] + Etot[j]) + hph * hph / (Ptot[j - 1] + Etot[j - 1]) / (hmh) + (2. - 1. / hdh) / (Ptot[j - 2] + Etot[j - 2]));
                }

                if ((RDIV - i) % 2 == 1) {
                    double fa = 1 / (Ptot[i] + Etot[i]);
                    double fb = 1 / (Ptot[i + 1] + Etot[i + 1]);
                    sum += ((Ptot[i] - Ptot[i + 1]) / 2.0) * (fa + fb);
                }

                hh = sum;
            }
            nu_gp[i] = nu_s - hh;
        }
        nu_gp[RDIV] = nu_s;
    }
}

/*C*/
/*************************************************************************/
void sphere(double s_gp[SDIV + 1], struct EOS *eosBM, struct stellar_properties *star_props, struct evolution_variables *ev, double *r_e_BM, double *s_e_BM, struct EOS *eosDM, struct stellar_properties *DM_props, double *r_e_DM, double *s_e_DM, double *massBM, double *massDM) {
    int s, m, n_nearest;
    bool is_DM;
    double r_is_s, r_is_final_big, r_is_final_small, r_final_big, r_final_small, m_final, lambda_s, nu_s, r_is_gp[RDIV + 1], lambda_gp[RDIV + 1], nu_gp[RDIV + 1], gama_mu_0[SDIV + 1], rho_mu_0[SDIV + 1], gama_eq, rho_eq, m_DM_final;

    double s_e = 0.5;

    /* The function TOV integrates the TOV equations. The function
           can be found in the file equil.c */
    is_DM = false;
    TOV(1, star_props, eosBM, r_is_gp, lambda_gp, nu_gp, &r_is_final_big, &r_final_big, &m_final, DM_props, eosDM, &r_is_final_small, &r_final_small, &m_DM_final, &is_DM);

    is_DM = false;
    TOV(2, star_props, eosBM, r_is_gp, lambda_gp, nu_gp, &r_is_final_big, &r_final_big, &m_final, DM_props, eosDM, &r_is_final_small, &r_final_small, &m_DM_final, &is_DM);

    is_DM = false;
    TOV(3, star_props, eosBM, r_is_gp, lambda_gp, nu_gp, &r_is_final_big, &r_final_big, &m_final, DM_props, eosDM, &r_is_final_small, &r_final_small, &m_DM_final, &is_DM);

    n_nearest = RDIV / 2;
    for (s = 1; s <= SDIV; s++) {
        r_is_s = r_is_final_big * (s_gp[s] / (1.0 - s_gp[s]));

        if (r_is_s < r_is_final_big) {
            lambda_s = interp(r_is_gp, lambda_gp, RDIV, r_is_s, &n_nearest);
            nu_s = interp(r_is_gp, nu_gp, RDIV, r_is_s, &n_nearest);
        } else {
            lambda_s = 2.0 * log(1.0 + (m_final + m_DM_final) / (2.0 * r_is_s));
            nu_s = log((1.0 - (m_final + m_DM_final) / (2.0 * r_is_s)) / (1.0 + (m_final + m_DM_final) / (2 * r_is_s)));
        }

        ev->gama[s][1] = nu_s + lambda_s;
        ev->rho[s][1] = nu_s - lambda_s;

        for (m = 1; m <= MDIV; m++) {
            ev->gama[s][m] = ev->gama[s][1];
            ev->rho[s][m] = ev->rho[s][1];
            ev->alpha[s][m] = (ev->gama[s][1] - ev->rho[s][1]) / 2.0;
            ev->omega[s][m] = 0.0;
        }

        gama_mu_0[s] = ev->gama[s][1]; /* gama at \mu=0 */
        rho_mu_0[s] = ev->rho[s][1];   /* rho at \mu=0 */
    }

    double r_e_big, r_e_small;
    n_nearest = SDIV / 2;
    gama_eq = interp(s_gp, gama_mu_0, SDIV, s_e, &n_nearest); /* gama at equator */
    rho_eq = interp(s_gp, rho_mu_0, SDIV, s_e, &n_nearest);   /* rho at equator */
    r_e_big = r_final_big * exp(0.5 * (rho_eq - gama_eq));

    double s_e_small = r_is_final_small / (r_is_final_small + r_is_final_big);
    gama_eq = interp(s_gp, gama_mu_0, SDIV, s_e_small, &n_nearest); /* gama at equator */
    rho_eq = interp(s_gp, rho_mu_0, SDIV, s_e_small, &n_nearest);   /* rho at equator */
    r_e_small = r_final_small * exp(0.5 * (rho_eq - gama_eq));

    *massBM = m_final * sqrt(KAPPA) * C * C / G / MSUN;
    *massDM = m_DM_final * sqrt(KAPPA) * C * C / G / MSUN;

    // If is_DM = true, the r_is_final_small is the DM radius.
    if (is_DM) {
        *r_e_BM = r_e_big;
        *r_e_DM = r_e_small;
        *s_e_BM = 0.5;
        *s_e_DM = s_e_small;
    } else {
        *r_e_BM = r_e_small;
        *r_e_DM = r_e_big;
        *s_e_BM = s_e_small;
        *s_e_DM = 0.5;
    }
}

/*C*/
/*************************************************************************/
/* Main iteration cycle for computation of the rotating star's metric    */
/*************************************************************************/

void log_divergence_info(double omega, double rho, double gama) {
    printf("Divergence detected: omega=%.2g, rho=%.2g, gama=%.2g\n", omega, rho, gama);
}

void compute_velocity_energy_pressure(double *s_gp, double *mu, struct evolution_variables *ev, struct stellar_properties *star_props, struct EOS *eos, struct DiffRotParams *DiffRotParams, double r_ratio, double r_ratio_other, double r_e, double s_e, double **Omega_diff, double **velocity_sq, double **enthalpy, double **pressure, double **energy, double gama_pole_h, double rho_pole_h, double *sin_theta, int verbose) {
    int n_nearest = eos->n_tab / 2;
    double sgp, mugp, enthalpy_value, common_term, rho0sm;

    for (int s = 1; s <= SDIV; s++) {
        sgp = s_gp[s];
        for (int m = 1; m <= MDIV; m++) {
            mugp = mu[m];

            // Compute velocity
            if ((r_ratio == 1.0 && r_ratio_other == 1.0) || s > SLIM) {
                velocity_sq[s][m] = 0.0;
            } else {
                velocity_sq[s][m] = SQ((Omega_diff[s][m] - ev->omega[s][m]) * (sgp / (1.0 - sgp)) * sin_theta[m] * exp(-ev->rho[s][m] * SQ(r_e)));
            }

            if (velocity_sq[s][m] >= 1.0) velocity_sq[s][m] = 0.0;

            // Compute enthalpy
            common_term = star_props->enthalpy_min + 0.5 * (SQ(r_e) * (gama_pole_h + rho_pole_h - ev->gama[s][m] - ev->rho[s][m]) - log(1.0 - velocity_sq[s][m]));

            double tmp;
            if (star_props->RotType == UNIFORM || r_ratio == 1) {
                enthalpy_value = common_term;
            } else if (star_props->RotType == DIFF) {
                tmp = -0.5 * SQ(DiffRotParams->A * (Omega_diff[s][m] - Omega_diff[1][1]));
                enthalpy_value = common_term + 0.5 * SQ(DiffRotParams->A * (Omega_diff[s][m] - Omega_diff[1][1]));
            } else if (star_props->RotType == U8) {
                double Omega_c = Omega_diff[1][1];

                double J = (Omega_diff[s][m] - ev->omega[s][m]) * SQ(sgp) * (1 - SQ(mugp)) * exp(-2 * SQ(r_e) * ev->rho[s][m]);
                J /= (SQ(1 - sgp) - (Omega_diff[s][m] - ev->omega[s][m]) * J);

                double tmp1 = SQ(DiffRotParams->A * Omega_c / 2) * (2 * SQ(DiffRotParams->A / DiffRotParams->B) * atan2(SQ(J), SQ(SQ(DiffRotParams->A) * Omega_c)) - sqrt(2.) * (atan2(SQ(DiffRotParams->A) * Omega_c - sqrt(2.) * J, SQ(DiffRotParams->A) * Omega_c) - atan2(SQ(DiffRotParams->A) * Omega_c + sqrt(2.) * J, SQ(DiffRotParams->A) * Omega_c)) + sqrt(2.) * atanh(sqrt(2.) * SQ(DiffRotParams->A) * Omega_c * J / (SQ(J) + SQ(SQ(DiffRotParams->A) * Omega_c))));

                tmp = -J * Omega_diff[s][m] + tmp1;
                enthalpy_value = common_term - J * Omega_diff[s][m] + tmp1;
            } else if (star_props->RotType == U9EXT) {
                double J = (Omega_diff[s][m] - ev->omega[s][m]) * SQ(sgp) * (1 - SQ(mugp)) * exp(-2 * SQ(r_e) * ev->rho[s][m]);
                J /= (SQ(1 - sgp) - (Omega_diff[s][m] - ev->omega[s][m]) * J);

                enthalpy_value = common_term - Omega_diff[1][1] * uryuExt_integral(J, Omega_diff[1][1], DiffRotParams->p, DiffRotParams->csi, DiffRotParams->A, DiffRotParams->B);
            }

            enthalpy[s][m] = enthalpy_value;

            // Compute energy and pressure
            if (enthalpy[s][m] <= star_props->enthalpy_min || sgp > s_e) {
                pressure[s][m] = 0.0;
                energy[s][m] = 0.0;
            } else {
                if (eos->type == TAB) {
                    pressure[s][m] = p_at_h(enthalpy_value, eos, &n_nearest);
                    energy[s][m] = e_at_p(pressure[s][m], eos, &n_nearest);
                } else {
                    rho0sm = pow(((eos->Gamma_P - 1.0) / eos->Gamma_P) * (exp(enthalpy_value) - 1.0), 1.0 / (eos->Gamma_P - 1.0));
                    pressure[s][m] = pow(rho0sm, eos->Gamma_P);
                    energy[s][m] = pressure[s][m] / (eos->Gamma_P - 1.0) + rho0sm;
                }
            }
        }  // m loop
    }  // s loop
}

void computeOmega(double s_gp[SDIV + 1], double mu[MDIV + 1], double **OmegaDiff, struct evolution_variables *ev, struct stellar_properties *props, struct DiffRotParams *DiffRotParam, double r_ratio, double *omega_mu_0, double *rho_mu_0, double s_e, double r_e, double gama_pole_h, double rho_pole_h, double gama_equator_h, double rho_equator_h, double *Omega_e, double *Omega_max, int n_of_it, int verbose, bool IsDM) {
    double Omega_c, Omega_h, omega_equator_h, term_in_Omega_h;
    double Jm_h, Je_h;
    double inf_j, sup_j, inf_Omega_e, sup_Omega_e, inf_Omega = 0, sup_Omega = 0;

    int n_nearest = SDIV / 2;

    int s, m;
    double const tol = DBL_EPSILON;

    int size = snprintf(NULL, 0, "Omega_diff_%s[%d][%d]", "BM", SDIV, MDIV) + 1;
    char *vName = malloc(size);
    if (vName == NULL) {
        fprintf(stderr, "Memory allocation failed for vName\n");
        return;
    }

    omega_equator_h = interp(s_gp, omega_mu_0, SDIV, s_e, &n_nearest);

    if (r_ratio == 1) {
        Omega_h = 0.0;
        omega_equator_h = 0.0;
        for (s = 1; s <= SDIV; s++) {
            for (m = 1; m <= MDIV; m++) {
                OmegaDiff[s][m] = Omega_h;
            }
        }
        *Omega_max = Omega_h;
    } else if (props->RotType == UNIFORM) {
        term_in_Omega_h = 1.0 - exp(SQ(r_e) * (gama_pole_h + rho_pole_h - gama_equator_h - rho_equator_h));
        if (term_in_Omega_h >= 0.0) {
            Omega_h = omega_equator_h + (1.0 - s_e) / s_e * exp(SQ(r_e) * rho_equator_h) * sqrt(term_in_Omega_h);
        } else {
            Omega_h = 0.0;
            printf("sqrt in Omega_h imaginary\n");
            exit(1);
        }

        // if (term_in_Omega_h < 0.){
        //     term_in_Omega_h = 0.5;
        //     DPRINT("SETTING v^2 < 0 to v^2 = 0.5! This should only happen in the first iterations, when the metric field change too much");
        // }
        // Omega_h = omega_equator_h + (1.0 - s_e) / s_e * exp(SQ(r_e) * rho_equator_h) * sqrt(term_in_Omega_h);

        for (s = 1; s <= SLIM; s++) {
            for (m = 1; m <= MDIV; m++) {
                OmegaDiff[s][m] = Omega_h;
            }
        }
        *Omega_max = Omega_h;
    } else if (props->RotType == DIFF) {
        find_boundary(call_diff_rotation, 0, 0.01, &inf_Omega_e, &sup_Omega_e, false, true, "Omega_e", 8, r_e, s_e, rho_equator_h, gama_equator_h, omega_equator_h, rho_pole_h, gama_pole_h, DiffRotParam->A);
        *Omega_e = brent_root_finder(call_diff_rotation, inf_Omega_e, sup_Omega_e, tol, "Omega_e", 8, r_e, s_e, rho_equator_h, gama_equator_h, omega_equator_h, rho_pole_h, gama_pole_h, DiffRotParam->A);

        // J_e = A^2 (Omega_c - Omega_e)
        Omega_c = (*Omega_e) + SQ(1.0 / DiffRotParam->A) * ((*Omega_e) - omega_equator_h) * SQ(s_e) * exp(-2.0 * SQ(r_e) * rho_equator_h) / (SQ(1.0 - s_e) - SQ(((*Omega_e) - omega_equator_h) * SQ(s_e) * exp(-SQ(r_e) * rho_equator_h)));

        Omega_h = omega_equator_h + exp(SQ(r_e) * rho_equator_h) * ((1.0 - s_e) / s_e) * sqrt(1.0 - exp(SQ(r_e) * (gama_pole_h + rho_pole_h - gama_equator_h - rho_equator_h)));

        s = 1;
        for (m = 1; m <= MDIV; m++)
            OmegaDiff[s][m] = Omega_c;

        m = MDIV;
        for (s = 1; s <= SLIM; s++)
            OmegaDiff[s][m] = Omega_c;

        for (s = 2; s <= SLIM; s++) {
            for (m = 1; m <= MDIV - 1; m++) {
                sprintf(vName, "Omega_diff_%s[%d][%d]", IsDM ? "DM" : "BM", s, m);

                int fail = 1;

                inf_Omega = OmegaDiff[s - 1][m] * 0.8;
                sup_Omega = OmegaDiff[s - 1][m] * 1.2;
                double llow = rotation_law(inf_Omega, r_e, ev->rho[s][m], ev->omega[s][m], s_gp[s], mu[m], Omega_c, DiffRotParam->A);
                double lhigh = rotation_law(sup_Omega, r_e, ev->rho[s][m], ev->omega[s][m], s_gp[s], mu[m], Omega_c, DiffRotParam->A);

                if (llow * lhigh < 0) {
                    fail = 0;
                }
                if (fail == 1) {
                    int fail = find_boundary(call_rotation_law, 0.0, 0.01, &inf_Omega, &sup_Omega, false, true, vName, 7, r_e, ev->rho[s][m], ev->omega[s][m], s_gp[s], mu[m], Omega_c, DiffRotParam->A);
                }
                OmegaDiff[s][m] = brent_root_finder(call_rotation_law, inf_Omega, sup_Omega, tol, vName, 7, r_e, ev->rho[s][m], ev->omega[s][m], s_gp[s], mu[m], Omega_c, DiffRotParam->A);
            }  // end of for m
        }  // end of for s

        if (verbose >= DEBUG) {
            FILE *fptro;
            if (IsDM)
                fptro = fopen("omegaDM.dat", "w");
            else
                fptro = fopen("omegaBM.dat", "w");
            m = 1;
            for (s = 1; s <= SLIM; s++) {
                fprintf(fptro, "%f %f %e\n", s_gp[s], mu[m], OmegaDiff[s][m] / r_e);
            }
            fclose(fptro);
        }

    } else if (props->RotType == U8) {
        double const p = 1;
        double const q = 3;

        // Compute Omega_e
        // if (n_of_it == 0) {
            find_boundary(call_uryu_diff_rotation, omega_equator_h + 10 * DBL_EPSILON, 0.01, &inf_Omega_e, &sup_Omega_e, false, true, "Omega_e", 10, r_e, s_e, rho_equator_h, gama_equator_h, omega_equator_h, rho_pole_h, gama_pole_h, DiffRotParam->A, DiffRotParam->B, DiffRotParam->lambda2);
        // } else {
        //     find_boundary(call_uryu_diff_rotation, *Omega_e * 1.05, -0.01, &sup_Omega_e, &inf_Omega_e, false, true, "Omega_e", 10, r_e, s_e, rho_equator_h, gama_equator_h, omega_equator_h, rho_pole_h, gama_pole_h, DiffRotParam->A, DiffRotParam->B, DiffRotParam->lambda2);
        // }

        // find_boundary(call_uryu_diff_rotation, n_of_it == 0 ? 0.5 : *Omega_e * 1.05, -0.01, &sup_Omega_e, &inf_Omega_e, false, true, "Omega_e", 10, r_e, s_e, rho_equator_h, gama_equator_h, omega_equator_h, rho_pole_h, gama_pole_h, DiffRotParam->A, DiffRotParam->B, DiffRotParam->lambda2);

        if (verbose >= DEBUG) {
            FILE *fptr;
            if (IsDM)
                fptr = fopen("omegaDM_e.dat", "w");
            else
                fptr = fopen("omegaBM_e.dat", "w");
            m = 1;
            double x = -0.5;
            while (x < 1) {
                double f1 = uryu_diff_rotation(x, r_e, s_e, rho_equator_h, gama_equator_h, omega_equator_h, rho_pole_h, gama_pole_h, DiffRotParam->A, DiffRotParam->B, DiffRotParam->lambda2);

                fprintf(fptr, "%f %e\n", x, f1);
                x += 0.01;
            }
            fclose(fptr);
        }

        // printf("Omega_e = %g\n", *Omega_e);

        *Omega_e = brent_root_finder(call_uryu_diff_rotation, inf_Omega_e, sup_Omega_e,
                                     tol, "Omega_e", 10, r_e, s_e, rho_equator_h, gama_equator_h, omega_equator_h, rho_pole_h, gama_pole_h, DiffRotParam->A, DiffRotParam->B, DiffRotParam->lambda2);
        Omega_c = *Omega_e / DiffRotParam->lambda2;

        // if (ev->omega[1][1] > Omega_c) {
        //     printf("Error! Angular velocity is smaller than frame dragging at s = 1, m = 1\n");
        //     printf("omega[1][1] = %g     Omega[1][1] = %g\n", ev->omega[1][1], Omega_c);
        // }
        // printf("Omega_e = %g\n", *Omega_e);
        // if (n_of_it == 0) exit(2);

        /* EQ A3 of  https://doi.org/10.1093/mnras/stab392 */
        Je_h = (*Omega_e - omega_equator_h) * SQ(s_e) * exp(-2 * SQ(r_e) * rho_equator_h) / (SQ(1 - s_e) - SQ(*Omega_e - omega_equator_h) * SQ(s_e) * exp(-2 * SQ(r_e) * rho_equator_h));

        // Compute Omega_max
        find_boundary(call_uryu_dOdj, DBL_EPSILON, 0.01, &inf_j, &sup_j, true, true, "J_max", 5, Omega_c, p, q, DiffRotParam->A, DiffRotParam->B);
        Jm_h = brent_root_finder(call_uryu_dOdj, inf_j, sup_j, tol, "J_max", 5, Omega_c, p, q, DiffRotParam->A, DiffRotParam->B);

        *Omega_max = Omega_c * (1 + pow(Jm_h / SQ(DiffRotParam->B) / Omega_c, p)) / (1 + pow(Jm_h / SQ(DiffRotParam->A) / Omega_c, p + q));

        Omega_h = omega_equator_h + exp(SQ(r_e) * rho_equator_h) * ((1 - s_e) / s_e) * sqrt(1.0 - exp(SQ(r_e) * (gama_pole_h + rho_pole_h - gama_equator_h - rho_equator_h)));

        // Compute Omega
        s = 1;
        for (m = 1; m <= MDIV; m++)
            OmegaDiff[s][m] = Omega_c;

        m = MDIV;
        for (s = 1; s <= SLIM; s++)
            OmegaDiff[s][m] = Omega_c;

        if (verbose >= DEBUG) {
            FILE *fptr;
            if (IsDM)
                fptr = fopen("omegaDM_diff_x.dat", "w");
            else
                fptr = fopen("omegaBM_diff_x.dat", "w");
            m = 1;
            for (s = 1; s <= SLIM; s++) {
                double x = -0.5;
                while (x < 1) {
                    double f1 = uryu_rotation_law(x, r_e, ev->rho[s][m], ev->omega[s][m], s_gp[s], mu[m], Omega_c, DiffRotParam->A, DiffRotParam->B);

                    fprintf(fptr, "%f  %e  %e  %d\n", s_gp[s], x, f1, m);
                    x += 0.01;
                }
            }
            fclose(fptr);
        }
        bool error = false;

        for (s = 2; s <= SLIM; s++) {
            for (m = 1; m <= MDIV - 1; m++) {
                // Create variable name for error handling
                sprintf(vName, "Omega_diff_%s[%d][%d]", IsDM ? "DM" : "BM", s, m);

                int fail = 1;

                // Check if broad range is enough to bound root
                // inf_Omega = OmegaDiff[s - 1][m] * 0.8;
                // sup_Omega = OmegaDiff[s - 1][m] * 1.2;
                // double llow = uryu_rotation_law(inf_Omega, r_e, ev->rho[s][m], ev->omega[s][m], s_gp[s], mu[m], Omega_c, DiffRotParam->A, DiffRotParam->B);
                // double lhigh = uryu_rotation_law(sup_Omega, r_e, ev->rho[s][m], ev->omega[s][m], s_gp[s], mu[m], Omega_c, DiffRotParam->A, DiffRotParam->B);

                // if (llow * lhigh < 0) {
                //     fail = 0;
                // }
                if (fail == 1) {
                    for (int i = 0; i < 2; i++) {
                        double step = 0.01 * pow(10, -(double)i);
                        fail = find_boundary(call_uryu_rotation_law, *Omega_max * 1.1, -step, &sup_Omega, &inf_Omega, false, false, vName, 8, r_e, ev->rho[s][m], ev->omega[s][m], s_gp[s], mu[m], Omega_c, DiffRotParam->A, DiffRotParam->B);
                        if (fail == 0) {
                            break;
                        }
                    }
                }

                if (fail == 1) {
                    // If there is an error, set it to the average to avoid crashing the code
                    // An error message will print at the end
                    error = true;
                    OmegaDiff[s][m] = (*Omega_max + *Omega_e) / 2;
                } else {
                    OmegaDiff[s][m] = brent_root_finder(call_uryu_rotation_law, inf_Omega, sup_Omega,
                                                        tol, vName, 8, r_e, ev->rho[s][m], ev->omega[s][m], s_gp[s], mu[m], Omega_c, DiffRotParam->A, DiffRotParam->B);
                    // OmegaDiff[s][m] = fmax(OmegaDiff[s][m], ev->omega[s][m]);
                }
            }  // end of for m
        }  // end of for s

        if (verbose >= DEBUG) {
            FILE *fptro;
            if (IsDM)
                fptro = fopen("omegaDM.dat", "w");
            else
                fptro = fopen("omegaBM.dat", "w");
            m = 1;
            for (s = 1; s <= SLIM; s++) {
                fprintf(fptro, "%f %f %e\n", s_gp[s], mu[m], OmegaDiff[s][m]);
            }
            fclose(fptro);
        }

    } else if (props->RotType == U9EXT) {
        if (verbose >= VERYVERBOSE) {
            DPRINT("\ns_e = %g\n", s_e);
            DPRINT("omega_equator_h = %g\n", omega_equator_h);
            DPRINT("Omega_e = %g\n", *Omega_e);
        }

        // 1) Find the new Omega_e from Eq. 25 of arXiv:2405.06609
        find_boundary(call_uryuExt_diff_rotation, omega_equator_h + 10 * DBL_EPSILON, 0.01, &inf_Omega_e, &sup_Omega_e, false, true, "Omega_e, computeOmega", 12, s_e, r_e, rho_equator_h, gama_equator_h, omega_equator_h, rho_pole_h, gama_pole_h, DiffRotParam->A, DiffRotParam->B, DiffRotParam->lambda2, DiffRotParam->p, DiffRotParam->csi);
        *Omega_e = brent_root_finder(call_uryuExt_diff_rotation, inf_Omega_e, sup_Omega_e,
                                     tol, "Omega_e, computeOmega", 12, s_e, r_e, rho_equator_h, gama_equator_h, omega_equator_h, rho_pole_h, gama_pole_h, DiffRotParam->A, DiffRotParam->B, DiffRotParam->lambda2, DiffRotParam->p, DiffRotParam->csi);
        if (verbose >= VERYVERBOSE) {
            DPRINT("Omega_e = %g inf_Omega_e = %g sup_Omega_e = %g\n", *Omega_e, inf_Omega_e, sup_Omega_e);
        }
        // 2) Update Je, Omega_c
        Je_h = (*Omega_e - omega_equator_h) * SQ(s_e) * exp(-2 * SQ(r_e) * rho_equator_h) / (SQ(1 - s_e) - SQ(*Omega_e - omega_equator_h) * SQ(s_e) * exp(-2. * SQ(r_e) * rho_equator_h));

        Omega_c = *Omega_e / DiffRotParam->lambda2;

        if (verbose >= VERYVERBOSE) {
            DPRINT("Je_h = %g\n", Je_h);

            DPRINT("Omega_c = %g\n", Omega_c);
        }
        // 3) Compute the angular momentum at maximum angular velocity from Eq. 49 of arXiv:2405.06609
        find_boundary(call_uryuExt_dOdj, DBL_EPSILON, 0.1, &inf_j, &sup_j, true, true, "Jm_h, computeOmega", 5, Omega_c, DiffRotParam->p, DiffRotParam->csi, DiffRotParam->A, DiffRotParam->B);
        Jm_h = brent_root_finder(call_uryuExt_dOdj, inf_j, sup_j, tol, "Jm_h, computeOmega", 5, Omega_c, DiffRotParam->p, DiffRotParam->csi, DiffRotParam->A, DiffRotParam->B);
        if (verbose >= VERYVERBOSE) {
            DPRINT("Jm_h = %g inf_j = %g sup_j = %g\n", Jm_h, inf_j, sup_j);
        }
        // Compute Omega_max from Eq. 48 of arXiv:2405.06609
        *Omega_max = Omega_c * (1 + pow(Jm_h / (SQ(DiffRotParam->B) * Omega_c), DiffRotParam->p)) * (1 - pow(Jm_h / (SQ(DiffRotParam->A) * Omega_c), DiffRotParam->csi));
        if (verbose >= VERYVERBOSE) {
            DPRINT("+ = %g    - = %g\n", Jm_h / (SQ(DiffRotParam->B) * Omega_c), Jm_h / (SQ(DiffRotParam->A) * Omega_c));
            DPRINT("Omega_max = %g\n", *Omega_max);
        }
        Omega_h = omega_equator_h + exp(SQ(r_e) * rho_equator_h) * ((1 - s_e) / s_e) * sqrt(1.0 - exp(SQ(r_e) * (gama_pole_h + rho_pole_h - gama_equator_h - rho_equator_h)));

        // 5) Fill Omega_diff solving Eq. 48 of arXiv:2405.06609 at each s and m
        s = 1;
        for (m = 1; m <= MDIV - 1; m++)
            OmegaDiff[s][m] = Omega_c;

        m = MDIV;
        for (s = 1; s <= SLIM; s++)
            OmegaDiff[s][m] = Omega_c;

        if (verbose >= DEBUG) {
            FILE *fptr;
            if (IsDM)
                fptr = fopen("omegaDM_diff_x.dat", "w");
            else
                fptr = fopen("omegaBM_diff_x.dat", "w");
            m = 1;
            for (s = 1; s <= SLIM; s++) {
                double x = -0.1;
                while (x < 1) {
                    double f1 = uryuExt_rotation_law(x, r_e,
                                                     ev->rho[s][m], ev->omega[s][m],
                                                     s_gp[s], mu[m], OmegaDiff[1][1],
                                                     DiffRotParam->A, DiffRotParam->B, DiffRotParam->p, DiffRotParam->csi);

                    fprintf(fptr, "%f  %e  %e  %d\n", s_gp[s], x, f1, m);
                    x += 0.001;
                }
            }
            fclose(fptr);
        }

        bool error = false;
        for (s = 2; s <= SLIM; s++) {
            for (m = 1; m <= MDIV - 1; m++) {
                // Create variable name for error handling
                sprintf(vName, "Omega_diff_%s[%d][%d]", IsDM ? "DM" : "BM", s, m);

                int fail = 1;

                // Check if broad range is enough to bound root
                inf_Omega = ev->omega[s][m] + 10 * DBL_EPSILON;
                sup_Omega = *Omega_max * 1.1;
                double llow = uryuExt_rotation_law(inf_Omega, r_e, ev->rho[s][m], ev->omega[s][m], s_gp[s], mu[m], Omega_c, DiffRotParam->A, DiffRotParam->B, DiffRotParam->p, DiffRotParam->csi);
                double lhigh = uryuExt_rotation_law(sup_Omega, r_e, ev->rho[s][m], ev->omega[s][m], s_gp[s], mu[m], Omega_c, DiffRotParam->A, DiffRotParam->B, DiffRotParam->p, DiffRotParam->csi);

                if (llow * lhigh < 0) {
                    fail = 0;
                }

                // If not, look for a good boundary
                if (fail == 1) {
                    find_boundary(call_uryuExt_rotation_law, *Omega_max * 1.1, -0.01, &sup_Omega, &inf_Omega, false, false, vName, 10, r_e, ev->rho[s][m], ev->omega[s][m], s_gp[s], mu[m], Omega_c, DiffRotParam->A, DiffRotParam->B, DiffRotParam->p, DiffRotParam->csi);
                }

                // Compute Omega_sm
                if (fail == 1) {
                    // If there is an error, set it to the average to avoid crashing the code
                    // An error message will print at the end
                    error = true;
                    OmegaDiff[s][m] = (*Omega_max + *Omega_e) / 2;
                } else {
                    OmegaDiff[s][m] = brent_root_finder(call_uryuExt_rotation_law, inf_Omega, sup_Omega,
                                                        tol, vName, 10, r_e, ev->rho[s][m], ev->omega[s][m], s_gp[s], mu[m], Omega_c, DiffRotParam->A, DiffRotParam->B, DiffRotParam->p, DiffRotParam->csi);
                }
            }  // end of for m
        }  // end of for s
        if (error) {
            if (IsDM) {
                printf("Error in computing DM!\n");
            } else {
                printf("Error in computing BM!\n");
            }
        }

        if (verbose >= DEBUG) {
            FILE *fptro;
            if (IsDM)
                fptro = fopen("omegaDM.dat", "w");
            else
                fptro = fopen("omegaBM.dat", "w");
            m = 1;
            for (s = 1; s <= SLIM; s++) {
                fprintf(fptro, "%f %f %e\n", s_gp[s], mu[m], OmegaDiff[s][m]);
            }
            fclose(fptro);
        }
    }

    // if (n_of_it == 4) exit(1);
    free(vName);
}

void updateAB(double s_gp[SDIV + 1], double mu[MDIV + 1], double **Omega_diff, struct evolution_variables *ev, struct stellar_properties *props, struct DiffRotParams *DiffRotParam, double r_ratio, double *omega_mu_0, double *rho_mu_0, double s_e, double r_e, double gama_pole_h, double rho_pole_h, double gama_equator_h, double rho_equator_h, int verbose, bool IsDM) {
    if (r_ratio == 1 || props->RotType == UNIFORM || props->RotType == DIFF) return;

    int n_nearest;
    double const tol = DBL_EPSILON;

    // Compute Je_h. This is agnostic to the rotation law, it only uses j = u_t u^\phi
    double Omega_equator[SLIM + 1], ss[SLIM + 1];
    for (int s = 1; s <= SLIM; s++) {
        ss[s] = s_gp[s];
        Omega_equator[s] = Omega_diff[s][1];
    }
    n_nearest = SDIV / 2;
    double omega_equator_h = interp(s_gp, omega_mu_0, SDIV, s_e, &n_nearest);
    n_nearest = SLIM / 2;
    double Omega_e = interp(ss, Omega_equator, SLIM, s_e, &n_nearest);
    double Je_h = (Omega_e - omega_equator_h) * SQ(s_e) * exp(-2 * SQ(r_e) * rho_equator_h) / (SQ(1 - s_e) - SQ(Omega_e - omega_equator_h) * SQ(s_e) * exp(-2 * SQ(r_e) * rho_equator_h));

    double Omega_c = Omega_e / DiffRotParam->lambda2;

    if (props->RotType == U8) {
        double p = 1;
        double q = 3;

        double inf_j, sup_j;
        find_boundary(call_uryu_dOdj, DBL_EPSILON, 0.01, &inf_j, &sup_j, true, true, "Jmax, Update_AB - U8", 5, Omega_c, p, q, DiffRotParam->A, DiffRotParam->B);
        double Jm_h = brent_root_finder(call_uryu_dOdj, inf_j, sup_j, tol, "Jmax, Update_AB - U8", 5, Omega_c, p, q, DiffRotParam->A, DiffRotParam->B);

        if (Je_h < pow(DiffRotParam->lambda1 / DiffRotParam->lambda2, 1.0 / q) * Jm_h) {
            printf("In %s:\n", IsDM ? "DM" : "BM");
            printf("s_e = %g   Omega_e = %g   Je_h = %g   Jm_h = %g\n", s_e, Omega_e, Je_h, Jm_h);
            printf("omega_equator_h = %g   rho_equator_h = %g\n", omega_equator_h, rho_equator_h);
            printf("Error diff rot: imaginary numbers for A2 and B2 \n");
            exit(1);
        }

        /* EQ A5 and A6 of  https://doi.org/10.1093/mnras/stab392 */
        DiffRotParam->A = sqrt((DiffRotParam->lambda2 / Omega_e) * pow(pow(Je_h * Jm_h, p) * (DiffRotParam->lambda2 * pow(Je_h, q) - DiffRotParam->lambda1 * pow(Jm_h, q)) / ((DiffRotParam->lambda1 - 1) * pow(Je_h, p) - (DiffRotParam->lambda2 - 1) * pow(Jm_h, p)), 1 / (p + q)));
        DiffRotParam->B = sqrt((Je_h * Jm_h * DiffRotParam->lambda2 / Omega_e) * pow((DiffRotParam->lambda2 * pow(Je_h, q) - DiffRotParam->lambda1 * pow(Jm_h, q)) / (DiffRotParam->lambda2 * (DiffRotParam->lambda1 - 1) * pow(Je_h, p + q) - DiffRotParam->lambda1 * (DiffRotParam->lambda2 - 1) * pow(Jm_h, p + q)), 1 / p));

    } else if (props->RotType == U9EXT) {
        // Find the parameters A, B
        // Solve analytically for B from Eq. 51 of arXiv:2405.06609 (error in sign of second to last term)

        double inf_j, sup_j;
        find_boundary(call_uryuExt_dOdj, DBL_EPSILON, 0.01, &inf_j, &sup_j, true, true, "Jmax, Update_AB - U9Ext", 5, Omega_c, DiffRotParam->p, DiffRotParam->csi, DiffRotParam->A, DiffRotParam->B);
        double Jm_h = brent_root_finder(call_uryuExt_dOdj, inf_j, sup_j, tol, "Jmax, Update_AB - U9Ext", 5, Omega_c, DiffRotParam->p, DiffRotParam->csi, DiffRotParam->A, DiffRotParam->B);

        find_boundary(call_uryuExt_Beq, DBL_EPSILON, 0.01, &inf_j, &sup_j, false, true, "B, Update_AB - U9Ext", 6, DiffRotParam->p, DiffRotParam->csi, DiffRotParam->lambda1, DiffRotParam->lambda2, Je_h, Jm_h);
        double x = brent_root_finder(call_uryuExt_Beq, inf_j, sup_j, tol, "B, Update_AB - U9Ext", 6, DiffRotParam->p, DiffRotParam->csi, DiffRotParam->lambda1, DiffRotParam->lambda2, Je_h, Jm_h);  // (B^2 Omega_c)^p
        DiffRotParam->B = pow(x, 1 / DiffRotParam->p) / Omega_c;

        // Solve for A from Eq. 52 of arXiv:2405.06609
        double num = pow(Jm_h, DiffRotParam->csi) * x + pow(Jm_h, DiffRotParam->p + DiffRotParam->csi);
        double den = x * (1. - DiffRotParam->lambda1) + pow(Jm_h, DiffRotParam->p);
        DiffRotParam->A = (1 / Omega_c) * pow(num / den, 1 / DiffRotParam->csi);

        if (DiffRotParam->A != DiffRotParam->A) {
            printf("In %s:\n", IsDM ? "DM" : "BM");
            printf("A is nan\n");
            printf("x = %f    num = %f    den = %f\n", x, num, den);
            printf("B = %f Je_h = %f Jm_h = %f Omega_c = %f \n\n", sqrt(DiffRotParam->B), Je_h, Jm_h, Omega_c);
            exit(EXIT_FAILURE);
        }

        // Compute A and B
        DiffRotParam->A = sqrt(DiffRotParam->A);
        DiffRotParam->B = sqrt(DiffRotParam->B);
    }
}

int spin(double s_gp[SDIV + 1], double mu[MDIV + 1], struct EOS *eosBM, struct stellar_properties *star_props, struct evolution_variables *ev, double accuracy, double cf, double r_ratio_BM, double *r_e_BM_new, double *Omega_BM, double s_e_BM, struct EOS *eosDM, struct stellar_properties *DM_props, double r_ratio_DM, double *r_e_DM_new, double *Omega_DM, double s_e_DM, int verbose, bool *outOfiter, bool counter, double *Omega_e_BM, struct DiffRotParams *DiffRotBM, double *Omega_e_DM, struct DiffRotParams *DiffRotDM, bool zero) {
    int m, s, n, k, n_of_it, i, j;
    double **D2_rho, **D2_gama, **D2_omega;
    float ***f_rho, ***f_gama;
    double sum_rho, sum_gama, sum_omega, r_e_BM_old, r_e_DM_old, d_gama_s, d_gama_m, d_rho_s, d_rho_m, d_omega_s, d_omega_m, d_gama_ss, d_gama_mm, d_gama_sm, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, m1, s1, s2, ea, rsm, gsm, omsm, esmBM, psmBM, v2smBM, OmegasmBM, OmegasmDM, esmDM, psmDM, v2smDM, mum, sgp, s_1, e_gsm, e_rsm, rho0sm, term_in_Omega_h, r_p_BM, s_p_BM, gama_pole_h_BM, gama_center_h, gama_equator_h_BM, rho_pole_h_BM, rho_center_h, rho_equator_h_BM, r_p_DM, s_p_DM, gama_pole_h_DM, gama_equator_h_DM, rho_pole_h_DM, rho_equator_h_DM, gama_mu_1[SDIV + 1], gama_mu_0[SDIV + 1], rho_mu_1[SDIV + 1], rho_mu_0[SDIV + 1], omega_mu_0[SDIV + 1], **da_dm, **dgds, **dgdm, **D1_rho, **D1_gama, **D1_omega, **S_gama, **S_rho, **S_omega, **f2n, **P_2n, **P1_2n_1, Omega_h_BM, Omega_h_DM, sin_theta[MDIV + 1], theta[MDIV + 1], r_e_BM, r_e_DM, Omega_max_BM, Omega_max_DM;
    double r_e, difBM = 1, difDM = 1;

    bool divergence_flag = false;

    n_of_it = 0;
    sum_rho = 0.0;
    sum_gama = 0.0;
    sum_omega = 0.0;

    /********************************************************************/
    /* PRECOMPUTE LEGENDRE POLYNOMIALS AND FUNCTIONS FOR INTEGRATION    */
    /********************************************************************/

    f2n = dmatrix(1, LMAX + 1, 1, SDIV);
    f_rho = f3tensor(1, SDIV, 1, LMAX + 1, 1, SDIV);
    f_gama = f3tensor(1, SDIV, 1, LMAX + 1, 1, SDIV);

    P_2n = dmatrix(1, MDIV, 1, LMAX + 1);
    P1_2n_1 = dmatrix(1, MDIV, 1, LMAX + 1);

    compute_polynomials(s_gp, mu, f_rho, f_gama, P_2n, P1_2n_1, sin_theta, theta);

    /********************************************************************/
    /* BEGINNING OF PROPER ALGORITHM                                    */
    /********************************************************************/
    S_gama = dmatrix(1, SDIV, 1, MDIV);
    S_rho = dmatrix(1, SDIV, 1, MDIV);
    S_omega = dmatrix(1, SDIV, 1, MDIV);
    D1_rho = dmatrix(1, LMAX + 1, 1, SDIV);
    D1_gama = dmatrix(1, LMAX + 1, 1, SDIV);
    D1_omega = dmatrix(1, LMAX + 1, 1, SDIV);
    D2_rho = dmatrix(1, SDIV, 1, LMAX + 1);
    D2_gama = dmatrix(1, SDIV, 1, LMAX + 1);
    D2_omega = dmatrix(1, SDIV, 1, LMAX + 1);
    da_dm = dmatrix(1, SDIV, 1, MDIV);
    dgds = dmatrix(1, SDIV, 1, MDIV);
    dgdm = dmatrix(1, SDIV, 1, MDIV);

    r_e_BM = *r_e_BM_new;
    r_e_DM = *r_e_DM_new;
    r_e = fmax(r_e_BM, r_e_DM);

    *Omega_e_BM = *Omega_e_BM * r_e;
    *Omega_e_DM = *Omega_e_DM * r_e;

    for (s = 1; s <= SDIV; s++) {
        for (m = 1; m <= MDIV; m++) {
            ev->Omega_diffBM[s][m] = ev->Omega_diffBM[s][m] * r_e;
            ev->Omega_diffDM[s][m] = ev->Omega_diffDM[s][m] * r_e;
        }
    }

    int over_counter = 0;
    double pdiffBM, pdiffDM;
    while ((difBM > accuracy || difDM > accuracy || over_counter < 4) || n_of_it < 2) {
        // Occasionally code is stuck in an infinite loop
        // Exit if it happens
        if (n_of_it > 100000) {
            *outOfiter = true;
            printf("Too many iterations!\n");
            break;
        }

        // Reduce convergence factor if too many iterations
        // are necessary. It usually means an error.
        if (r_ratio_BM == 1 && r_ratio_DM == 1) {
            cf = 1;
        } else if (n_of_it < 500) {
            cf = 0.5;
        } else if (n_of_it < 600) {
            cf = 0.1;
        } else if (n_of_it < 800) {
            cf = 0.01;
        } else if (n_of_it < 1000) {
            cf = 0.001;
        } else {
            cf = 0.0001;
        }

        // Rescale potentials and construct arrays with
        // the potentials along the equatorial and polar directions
        for (s = 1; s <= SDIV; s++) {
            for (m = 1; m <= MDIV; m++) {
                ev->rho[s][m] /= SQ(r_e);
                ev->gama[s][m] /= SQ(r_e);
                ev->alpha[s][m] /= SQ(r_e);
                ev->omega[s][m] *= r_e;
            }
            rho_mu_0[s] = ev->rho[s][1];
            gama_mu_0[s] = ev->gama[s][1];
            omega_mu_0[s] = ev->omega[s][1];
            rho_mu_1[s] = ev->rho[s][MDIV];
            gama_mu_1[s] = ev->gama[s][MDIV];
        }

        gama_center_h = ev->gama[1][1];
        rho_center_h = ev->rho[1][1];

        // BM
        // Compute new r_e.
        r_e_BM_old = r_e_BM;
        r_p_BM = r_ratio_BM * r_e_BM;
        s_e_BM = r_e_BM / (r_e_BM + r_e);
        s_p_BM = r_p_BM / (r_p_BM + r_e);

        int n_nearest = SDIV / 2;
        gama_pole_h_BM = interp(s_gp, gama_mu_1, SDIV, s_p_BM, &n_nearest);
        gama_equator_h_BM = interp(s_gp, gama_mu_0, SDIV, s_e_BM, &n_nearest);

        rho_pole_h_BM = interp(s_gp, rho_mu_1, SDIV, s_p_BM, &n_nearest);
        rho_equator_h_BM = interp(s_gp, rho_mu_0, SDIV, s_e_BM, &n_nearest);

        if (r_ratio_BM != 1 & n_of_it > 0) {
            updateAB(s_gp, mu, ev->Omega_diffBM, ev, star_props, DiffRotBM, r_ratio_BM, omega_mu_0, rho_mu_0, s_e_BM, r_e, gama_pole_h_BM, rho_pole_h_BM, gama_equator_h_BM, rho_equator_h_BM, verbose, false);
        }

        double tmp = 2.0 * star_props->h_center / ((gama_pole_h_BM + rho_pole_h_BM - gama_center_h - rho_center_h) * SQ(r_e));
        r_e_BM = r_e_BM * sqrt(tmp);

        // DM
        // Compute new r_e.
        if (!zero) {
            r_e_DM_old = r_e_DM;
            r_p_DM = r_ratio_DM * r_e_DM;
            s_e_DM = r_e_DM / (r_e_DM + r_e);
            s_p_DM = r_p_DM / (r_p_DM + r_e);

            n_nearest = SDIV / 2;
            gama_pole_h_DM = interp(s_gp, gama_mu_1, SDIV, s_p_DM, &n_nearest);
            gama_equator_h_DM = interp(s_gp, gama_mu_0, SDIV, s_e_DM, &n_nearest);

            rho_pole_h_DM = interp(s_gp, rho_mu_1, SDIV, s_p_DM, &n_nearest);
            rho_equator_h_DM = interp(s_gp, rho_mu_0, SDIV, s_e_DM, &n_nearest);

            if (n_of_it > 0)
                updateAB(s_gp, mu, ev->Omega_diffDM, ev, DM_props, DiffRotDM, r_ratio_DM, omega_mu_0, rho_mu_0, s_e_DM, r_e, gama_pole_h_DM, rho_pole_h_DM, gama_equator_h_DM, rho_equator_h_DM, verbose, true);

            double tmp = 2.0 * (DM_props->h_center - DM_props->enthalpy_min) / ((gama_pole_h_DM + rho_pole_h_DM - gama_center_h - rho_center_h) * SQ(r_e));
            r_e_DM = r_e_DM * sqrt(tmp);
        }

        // Update r_e
        r_e = fmax(r_e_BM, r_e_DM);

        // BM
        r_p_BM = r_ratio_BM * r_e_BM;
        s_e_BM = r_e_BM / (r_e_BM + r_e);
        s_p_BM = r_p_BM / (r_p_BM + r_e);

        gama_pole_h_BM = interp(s_gp, gama_mu_1, SDIV, s_p_BM, &n_nearest);
        gama_equator_h_BM = interp(s_gp, gama_mu_0, SDIV, s_e_BM, &n_nearest);

        rho_pole_h_BM = interp(s_gp, rho_mu_1, SDIV, s_p_BM, &n_nearest);
        rho_equator_h_BM = interp(s_gp, rho_mu_0, SDIV, s_e_BM, &n_nearest);

        // Compute angular velocity.
        computeOmega(s_gp, mu, ev->Omega_diffBM, ev, star_props, DiffRotBM, r_ratio_BM, omega_mu_0, rho_mu_0, s_e_BM, r_e, gama_pole_h_BM, rho_pole_h_BM, gama_equator_h_BM, rho_equator_h_BM, Omega_e_BM, &Omega_max_BM, n_of_it, verbose, false);

        // Compute velocity, energy density and pressure.
        compute_velocity_energy_pressure(s_gp, mu, ev, star_props, eosBM, DiffRotBM, r_ratio_BM, r_ratio_DM, r_e, s_e_BM, ev->Omega_diffBM, ev->velocity_sqBM, ev->enthalpyBM, ev->pressureBM, ev->energyBM, gama_pole_h_BM, rho_pole_h_BM, sin_theta, verbose);

        // DM
        if (!zero) {
            r_p_DM = r_ratio_DM * r_e_DM;
            s_e_DM = r_e_DM / (r_e_DM + r_e);
            s_p_DM = r_p_DM / (r_p_DM + r_e);

            gama_pole_h_DM = interp(s_gp, gama_mu_1, SDIV, s_p_DM, &n_nearest);
            gama_equator_h_DM = interp(s_gp, gama_mu_0, SDIV, s_e_DM, &n_nearest);

            rho_pole_h_DM = interp(s_gp, rho_mu_1, SDIV, s_p_DM, &n_nearest);
            rho_equator_h_DM = interp(s_gp, rho_mu_0, SDIV, s_e_DM, &n_nearest);

            // Compute angular velocity.
            computeOmega(s_gp, mu, ev->Omega_diffDM, ev, DM_props, DiffRotDM, r_ratio_DM, omega_mu_0, rho_mu_0, s_e_DM, r_e, gama_pole_h_DM, rho_pole_h_DM, gama_equator_h_DM, rho_equator_h_DM, Omega_e_DM, &Omega_max_DM, n_of_it, verbose, true);

            // Compute velocity, energy density and pressure.
            compute_velocity_energy_pressure(s_gp, mu, ev, DM_props, eosDM, DiffRotDM, r_ratio_DM, r_ratio_BM, r_e, s_e_DM, ev->Omega_diffDM, ev->velocity_sqDM, ev->enthalpyDM, ev->pressureDM, ev->energyDM, gama_pole_h_DM, rho_pole_h_DM, sin_theta, verbose);
        }

        // Rescale back metric potentials (except omega)
        for (s = 1; s <= SDIV; s++) {
            for (m = 1; m <= MDIV; m++) {
                ev->rho[s][m] *= SQ(r_e);
                ev->gama[s][m] *= SQ(r_e);
                ev->alpha[s][m] *= SQ(r_e);
            }
        }

        // Compute metric source potentials
        for (s = 1; s <= SDIV; s++) {
            for (m = 1; m <= MDIV; m++) {
                rsm = ev->rho[s][m];
                gsm = ev->gama[s][m];
                omsm = ev->omega[s][m];
                esmBM = ev->energyBM[s][m];
                psmBM = ev->pressureBM[s][m];
                esmDM = ev->energyDM[s][m];
                psmDM = ev->pressureDM[s][m];
                OmegasmBM = ev->Omega_diffBM[s][m];
                OmegasmDM = ev->Omega_diffDM[s][m];
                e_gsm = exp(0.5 * gsm);
                e_rsm = exp(-rsm);
                v2smBM = ev->velocity_sqBM[s][m];
                v2smDM = ev->velocity_sqDM[s][m];
                mum = mu[m];
                m1 = 1.0 - SQ(mum);
                sgp = s_gp[s];
                s_1 = 1.0 - sgp;
                s1 = sgp * s_1;
                s2 = SQ(sgp / s_1);

                ea = 16.0 * PI * exp(2.0 * ev->alpha[s][m]) * SQ(r_e);

                if (s == 1) {
                    d_gama_s = 0.0;
                    d_gama_m = 0.0;
                    d_rho_s = 0.0;
                    d_rho_m = 0.0;
                    d_omega_s = 0.0;
                    d_omega_m = 0.0;
                } else {
                    d_gama_s = deriv_s(ev->gama, s, m);
                    d_gama_m = deriv_m(ev->gama, s, m);
                    d_rho_s = deriv_s(ev->rho, s, m);
                    d_rho_m = deriv_m(ev->rho, s, m);
                    d_omega_s = deriv_s(ev->omega, s, m);
                    d_omega_m = deriv_m(ev->omega, s, m);
                }

                S_rho[s][m] = e_gsm * (0.5 * ea * (esmBM + psmBM) * s2 * (1.0 + v2smBM) / (1.0 - v2smBM) + 0.5 * ea * (esmDM + psmDM) * s2 * (1.0 + v2smDM) / (1.0 - v2smDM) + s2 * m1 * SQ(e_rsm) * (SQ(s1 * d_omega_s) + m1 * SQ(d_omega_m)) + s1 * d_gama_s - mum * d_gama_m + 0.5 * rsm * (ea * (psmBM + psmDM) * s2 - s1 * d_gama_s * (0.5 * s1 * d_gama_s + 1.0) - d_gama_m * (0.5 * m1 * d_gama_m - mum)));

                S_gama[s][m] = e_gsm * (ea * (psmBM + psmDM) * s2 + 0.5 * gsm * (ea * (psmBM + psmDM) * s2 - 0.5 * SQ(s1 * d_gama_s) - 0.5 * m1 * SQ(d_gama_m)));

                S_omega[s][m] = e_gsm * e_rsm * (-ea * (OmegasmBM - omsm) * (esmBM + psmBM) * s2 / (1.0 - v2smBM) - ea * (OmegasmDM - omsm) * (esmDM + psmDM) * s2 / (1.0 - v2smDM) + omsm * (-0.5 * ea * (((1.0 + v2smBM) * esmBM + 2.0 * v2smBM * psmBM) / (1.0 - v2smBM)) * s2 - 0.5 * ea * (((1.0 + v2smDM) * esmDM + 2.0 * v2smDM * psmDM) / (1.0 - v2smDM)) * s2 - s1 * (2 * d_rho_s + 0.5 * d_gama_s) + mum * (2 * d_rho_m + 0.5 * d_gama_m) + 0.25 * SQ(s1) * (4 * SQ(d_rho_s) - SQ(d_gama_s)) + 0.25 * m1 * (4 * SQ(d_rho_m) - SQ(d_gama_m)) - m1 * SQ(e_rsm) * (SQ(SQ(sgp) * d_omega_s) + s2 * m1 * SQ(d_omega_m))));
            }
        }

        /* ANGULAR INTEGRATION */

        n = 0;
        for (k = 1; k <= SDIV; k++) {
            for (m = 1; m <= MDIV - 2; m += 2) {
                sum_rho += (DM / 3.0) * (P_2n[m][n + 1] * S_rho[k][m] + 4.0 * P_2n[m + 1][n + 1] * S_rho[k][m + 1] + P_2n[m + 2][n + 1] * S_rho[k][m + 2]);
            }
            D1_rho[n + 1][k] = sum_rho;
            D1_gama[n + 1][k] = 0.0;
            D1_omega[n + 1][k] = 0.0;
            sum_rho = 0.0;
        }

        for (n = 1; n <= LMAX; n++)
            for (k = 1; k <= SDIV; k++) {
                for (m = 1; m <= MDIV - 2; m += 2) {
                    sum_rho += (DM / 3.0) * (P_2n[m][n + 1] * S_rho[k][m] + 4.0 * P_2n[m + 1][n + 1] * S_rho[k][m + 1] + P_2n[m + 2][n + 1] * S_rho[k][m + 2]);
                    sum_gama += (DM / 3.0) * (sin((2.0 * n - 1.0) * theta[m]) * S_gama[k][m] + 4.0 * sin((2.0 * n - 1.0) * theta[m + 1]) * S_gama[k][m + 1] + sin((2.0 * n - 1.0) * theta[m + 2]) * S_gama[k][m + 2]);
                    sum_omega += (DM / 3.0) * (sin_theta[m] * P1_2n_1[m][n + 1] * S_omega[k][m] + 4.0 * sin_theta[m + 1] * P1_2n_1[m + 1][n + 1] * S_omega[k][m + 1] + sin_theta[m + 2] * P1_2n_1[m + 2][n + 1] * S_omega[k][m + 2]);
                }

                D1_rho[n + 1][k] = sum_rho;
                D1_gama[n + 1][k] = sum_gama;
                D1_omega[n + 1][k] = sum_omega;
                sum_rho = 0.0;
                sum_gama = 0.0;
                sum_omega = 0.0;
            }

        /* RADIAL INTEGRATION */

        n = 0;
        for (s = 1; s <= SDIV; s++) {
            for (k = 1; k <= SDIV - 2; k += 2) {
                sum_rho += (DS / 3.0) * (f_rho[s][n + 1][k] * D1_rho[n + 1][k] + 4.0 * f_rho[s][n + 1][k + 1] * D1_rho[n + 1][k + 1] + f_rho[s][n + 1][k + 2] * D1_rho[n + 1][k + 2]);
            }
            D2_rho[s][n + 1] = sum_rho;
            D2_gama[s][n + 1] = 0.0;
            D2_omega[s][n + 1] = 0.0;
            sum_rho = 0.0;
        }

        for (s = 1; s <= SDIV; s++) {
            for (n = 1; n <= LMAX; n++) {
                for (k = 1; k <= SDIV - 2; k += 2) {
                    sum_rho += (DS / 3.0) * (f_rho[s][n + 1][k] * D1_rho[n + 1][k] + 4.0 * f_rho[s][n + 1][k + 1] * D1_rho[n + 1][k + 1] + f_rho[s][n + 1][k + 2] * D1_rho[n + 1][k + 2]);
                    sum_gama += (DS / 3.0) * (f_gama[s][n + 1][k] * D1_gama[n + 1][k] + 4.0 * f_gama[s][n + 1][k + 1] * D1_gama[n + 1][k + 1] + f_gama[s][n + 1][k + 2] * D1_gama[n + 1][k + 2]);

                    if (k < s && k + 2 <= s) {
                        sum_omega += (DS / 3.0) * (f_rho[s][n + 1][k] * D1_omega[n + 1][k] + 4.0 * f_rho[s][n + 1][k + 1] * D1_omega[n + 1][k + 1] + f_rho[s][n + 1][k + 2] * D1_omega[n + 1][k + 2]);
                    } else {
                        if (k >= s) {
                            sum_omega += (DS / 3.0) * (f_gama[s][n + 1][k] * D1_omega[n + 1][k] + 4.0 * f_gama[s][n + 1][k + 1] * D1_omega[n + 1][k + 1] + f_gama[s][n + 1][k + 2] * D1_omega[n + 1][k + 2]);
                        } else {
                            sum_omega += (DS / 3.0) * (f_rho[s][n + 1][k] * D1_omega[n + 1][k] + 4.0 * f_rho[s][n + 1][k + 1] * D1_omega[n + 1][k + 1] + f_gama[s][n + 1][k + 2] * D1_omega[n + 1][k + 2]);  // Why last term is is f_gama?
                        }
                    }
                }
                D2_rho[s][n + 1] = sum_rho;
                D2_gama[s][n + 1] = sum_gama;
                D2_omega[s][n + 1] = sum_omega;
                sum_rho = 0.0;
                sum_gama = 0.0;
                sum_omega = 0.0;
            }
        }

        /* SUMMATION OF COEFFICIENTS */
        for (s = 1; s <= SDIV; s++) {
            for (m = 1; m <= MDIV; m++) {
                gsm = ev->gama[s][m];
                rsm = ev->rho[s][m];
                omsm = ev->omega[s][m];
                e_gsm = exp(-0.5 * gsm);
                e_rsm = exp(rsm);
                temp1 = sin_theta[m];
                sum_rho += -e_gsm * P_2n[m][0 + 1] * D2_rho[s][0 + 1];

                for (n = 1; n <= LMAX; n++) {
                    sum_rho += -e_gsm * P_2n[m][n + 1] * D2_rho[s][n + 1];

                    if (m == MDIV) {
                        sum_omega += 0.5 * e_rsm * e_gsm * D2_omega[s][n + 1];
                        sum_gama += -(2.0 / PI) * e_gsm * D2_gama[s][n + 1];
                    } else {
                        sum_omega += -e_rsm * e_gsm * (P1_2n_1[m][n + 1] / (2.0 * n * (2.0 * n - 1.0) * temp1)) * D2_omega[s][n + 1];
                        sum_gama += -(2.0 / PI) * e_gsm * (sin((2.0 * n - 1.0) * theta[m]) / ((2.0 * n - 1.0) * temp1)) * D2_gama[s][n + 1];
                    }
                }

                ev->rho[s][m] = rsm + cf * (sum_rho - rsm);
                ev->gama[s][m] = gsm + cf * (sum_gama - gsm);
                ev->omega[s][m] = omsm + cf * (sum_omega - omsm);

                sum_omega = 0.0;
                sum_rho = 0.0;
                sum_gama = 0.0;
            }
        }

        /* CHECK FOR DIVERGENCE */

        if (fabs(ev->omega[2][1]) > 100.0 || fabs(ev->rho[2][1]) > 100.0 || fabs(ev->gama[2][1]) > 300.0) {
            log_divergence_info(ev->omega[2][1], ev->rho[2][1], ev->gama[2][1]);
            divergence_flag = true;
            break;
        }

        /* TREAT SPHERICAL CASE */

        if ((r_ratio_BM == 1.0) && (r_ratio_DM == 1.0)) {
            for (s = 1; s <= SDIV; s++)
                for (m = 1; m <= MDIV; m++) {
                    ev->rho[s][m] = ev->rho[s][1];
                    ev->gama[s][m] = ev->gama[s][1];
                    ev->omega[s][m] = 0.0;
                }
        }

        /* TREAT INFINITY WHEN SMAX=1.0 */
        if (SMAX == 1.0) {
            for (m = 1; m <= MDIV; m++) {
                ev->rho[SDIV][m] = 0.0;
                ev->gama[SDIV][m] = 0.0;
                ev->omega[SDIV][m] = 0.0;
            }
        }

        /* COMPUTE FIRST ORDER DERIVATIVES OF GAMA */

        for (s = 1; s <= SDIV; s++) {
            for (m = 1; m <= MDIV; m++) {
                dgds[s][m] = deriv_s(ev->gama, s, m);
                dgdm[s][m] = deriv_m(ev->gama, s, m);
            }
        }

        /* ALPHA */

        if ((r_ratio_BM == 1.0) && (r_ratio_DM == 1.0)) {
            for (s = 1; s <= SDIV; s++)
                for (m = 1; m <= MDIV; m++)
                    da_dm[s][m] = 0.0;
        } else {
            for (s = 2; s <= SDIV; s++) {
                for (m = 1; m <= MDIV; m++) {
                    if (s == 2) {
                        da_dm[1][m] = 0.0;
                    }
                    sgp = s_gp[s];
                    s1 = sgp * (1.0 - sgp);
                    mum = mu[m];
                    m1 = 1.0 - SQ(mum);

                    d_gama_s = dgds[s][m];
                    d_gama_m = dgdm[s][m];
                    d_rho_s = deriv_s(ev->rho, s, m);
                    d_rho_m = deriv_m(ev->rho, s, m);
                    d_omega_s = deriv_s(ev->omega, s, m);
                    d_omega_m = deriv_m(ev->omega, s, m);
                    d_gama_ss = s1 * deriv_s(dgds, s, m) + (1.0 - 2.0 * sgp) * d_gama_s;
                    d_gama_mm = m1 * deriv_m(dgdm, s, m) - 2.0 * mum * d_gama_m;
                    d_gama_sm = deriv_sm(ev->gama, s, m);

                    temp1 = 2.0 * SQ(sgp) * (sgp / (1.0 - sgp)) * m1 * d_omega_s * d_omega_m * (1.0 + s1 * d_gama_s) - (SQ(SQ(sgp) * d_omega_s) - SQ(sgp * d_omega_m / (1.0 - sgp)) * m1) * (-mum + m1 * d_gama_m);

                    temp2 = 1.0 / (m1 * SQ(1.0 + s1 * d_gama_s) + SQ(-mum + m1 * d_gama_m));

                    temp3 = s1 * d_gama_ss + SQ(s1 * d_gama_s);

                    temp4 = d_gama_m * (-mum + m1 * d_gama_m);

                    temp5 = (SQ(s1 * (d_rho_s + d_gama_s)) - m1 * SQ(d_rho_m + d_gama_m)) * (-mum + m1 * d_gama_m);

                    temp6 = s1 * m1 * (0.5 * (d_rho_s + d_gama_s) * (d_rho_m + d_gama_m) + d_gama_sm + d_gama_s * d_gama_m) * (1.0 + s1 * d_gama_s);

                    temp7 = s1 * mum * d_gama_s * (1.0 + s1 * d_gama_s);

                    temp8 = m1 * exp(-2 * ev->rho[s][m]);

                    da_dm[s][m] = -0.5 * (d_rho_m + d_gama_m) - temp2 * (0.5 * (temp3 - d_gama_mm - temp4) * (-mum + m1 * d_gama_m) + 0.25 * temp5 - temp6 + temp7 + 0.25 * temp8 * temp1);
                }
            }
        }

        for (s = 1; s <= SDIV; s++) {
            ev->alpha[s][1] = 0.0;
            for (m = 1; m <= MDIV - 1; m++) {
                ev->alpha[s][m + 1] = ev->alpha[s][m] + 0.5 * DM * (da_dm[s][m + 1] + da_dm[s][m]);
            }
        }

        for (s = 1; s <= SDIV; s++) {
            for (m = 1; m <= MDIV; m++) {
                ev->alpha[s][m] += -ev->alpha[s][MDIV] + 0.5 * (ev->gama[s][MDIV] - ev->rho[s][MDIV]);

                if (ev->alpha[s][m] >= 300.0) {
                    printf("alpha[%d][%d] = %g\n", s, m, ev->alpha[s][m]);
                    divergence_flag = true;
                    break;
                }
                ev->omega[s][m] /= r_e;
            }
            if (divergence_flag) {
                printf("Breaking loop\n");
                break;
            }
        }

        if (divergence_flag) {
            printf("Stopping!\n");
            break;
        }

        if (SMAX == 1.0) {
            for (m = 1; m <= MDIV; m++)
                ev->alpha[SDIV][m] = 0.0;
        }

        pdiffBM = difBM;
        pdiffDM = difDM;
        difBM = fabs(r_e_BM_old - r_e_BM) / r_e_BM;
        difDM = zero ? 0 : fabs(r_e_DM_old - r_e_DM) / r_e_DM;

        if (difBM < accuracy && difDM < accuracy) {
            over_counter++;
        } else {
            over_counter = 0;
        }

        if (verbose >= VERBOSE) {
            if (zero) {
                double l1 = Omega_max_BM / (DBL_EPSILON + ev->Omega_diffBM[1][1]);
                printf("it = %d  cf = %g  r_e_BM = %f  difBM = %4.3e lambda_1_BM = %g \n", n_of_it, cf, r_e_BM, difBM, l1);
            } else {
                double l1BM = Omega_max_BM / (DBL_EPSILON + ev->Omega_diffBM[1][1]);
                double l1DM = Omega_max_DM / (DBL_EPSILON + ev->Omega_diffDM[1][1]);
                printf("it = %d  cf = %g  r_e_BM = %f  r_e_DM = %f  difBM = %4.3e  difDM = %4.3e  lambda_1_BM = %g  lambda_1_DM = %g\n", n_of_it, cf, r_e_BM, r_e_DM, difBM, difDM, l1BM, l1DM);
            }
        }

        // char output_path[1024], input_filename[1024];
        // snprintf(output_path, sizeof(output_path), "%s/%s", ".", "test_nonSmooth.h5");
        // snprintf(input_filename, sizeof(input_filename), "%s", "/home/lorenzo/phd/Bifluid/RNS/rnsplusplus_diff/test/wrongProf/config.d");

        // for (s = 1; s <= SLIM; s++) {
        //     for (m = 1; m <= MDIV; m++) {
        //         ev->Omega_diffBM[s][m] = ev->Omega_diffBM[s][m] / r_e;
        //         ev->Omega_diffDM[s][m] = ev->Omega_diffDM[s][m] / r_e;
        //     }
        // }

        // hdf5_save_var(s_gp, SDIV, mu, MDIV, output_path, ev, star_props,
        //               DM_props, eosBM, eosDM, &r_ratio_BM, &r_ratio_DM,
        //               Omega_BM, Omega_DM, &r_e_BM, &r_e_DM,
        //               DiffRotBM, DiffRotDM, input_filename);

        // for (s = 1; s <= SLIM; s++) {
        //     for (m = 1; m <= MDIV; m++) {
        //         ev->Omega_diffBM[s][m] = ev->Omega_diffBM[s][m] * r_e;
        //         ev->Omega_diffDM[s][m] = ev->Omega_diffDM[s][m] * r_e;
        //     }
        // }

        // if (n_of_it == 4) exit(2);

        n_of_it++;
    } /* end while */

    free_dmatrix(S_gama, 1, SDIV, 1, MDIV);
    free_dmatrix(S_rho, 1, SDIV, 1, MDIV);
    free_dmatrix(S_omega, 1, SDIV, 1, MDIV);
    free_dmatrix(D1_rho, 1, LMAX + 1, 1, SDIV);
    free_dmatrix(D1_gama, 1, LMAX + 1, 1, SDIV);
    free_dmatrix(D1_omega, 1, LMAX + 1, 1, SDIV);
    free_dmatrix(D2_rho, 1, SDIV, 1, LMAX + 1);
    free_dmatrix(D2_gama, 1, SDIV, 1, LMAX + 1);
    free_dmatrix(D2_omega, 1, SDIV, 1, LMAX + 1);
    free_dmatrix(da_dm, 1, SDIV, 1, MDIV);
    free_dmatrix(dgds, 1, SDIV, 1, MDIV);
    free_dmatrix(dgdm, 1, SDIV, 1, MDIV);
    free_f3tensor(f_rho, 1, SDIV, 1, LMAX + 1, 1, SDIV);
    free_f3tensor(f_gama, 1, SDIV, 1, LMAX + 1, 1, SDIV);
    free_dmatrix(f2n, 1, LMAX + 1, 1, SDIV);
    free_dmatrix(P_2n, 1, MDIV, 1, LMAX + 1);
    free_dmatrix(P1_2n_1, 1, MDIV, 1, LMAX + 1);

    if (divergence_flag) {
        exit(EXIT_FAILURE);
    }

    /* COMPUTE OMEGA */

    (*Omega_BM) = ev->Omega_diffBM[1][1] / r_e;
    if (eosBM->type == TAB)
        (*Omega_BM) *= C / sqrt(KAPPA);

    (*Omega_DM) = ev->Omega_diffDM[1][1] / r_e;
    if (eosDM->type == TAB)
        (*Omega_DM) *= C / sqrt(KAPPA);

    /* Rescale back Omega_diff[s][m] */
    for (s = 1; s <= SLIM; s++) {
        for (m = 1; m <= MDIV; m++) {
            ev->Omega_diffBM[s][m] = ev->Omega_diffBM[s][m] / r_e;
            ev->Omega_diffDM[s][m] = ev->Omega_diffDM[s][m] / r_e;
        }
    }

    *Omega_e_BM = *Omega_e_BM / r_e;
    *Omega_e_DM = *Omega_e_DM / r_e;

    /* UPDATE r_e_new */

    *r_e_BM_new = r_e_BM;
    *r_e_DM_new = r_e_DM;

    return 0;
}