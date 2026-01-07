#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "consts.h"
#include "equil.h"
#include "equil_util.h"
#include "nrutil.h"
#include "output.h"

void set_default_values(double *cf, double *accuracy, double *r_ratio_BM, double *r_ratio_DM, double *fdm_target, double *MDM, double *fdm_tol, int *verbose, struct stellar_properties *star_props, struct stellar_properties *DM_props, struct EOS *eosBM, struct EOS *eosDM, struct DiffRotParams *DiffRotBM, struct DiffRotParams *DiffRotDM, char *input_filename, char *output2d_filename, char *outputS_filename, char *output_name, char *id_file) {
    *cf = 0.8;
    *accuracy = 1e-12;
    *r_ratio_BM = 1;
    *r_ratio_DM = 1;
    *fdm_target = 0.05;

    *MDM = 1000;  // MeV - DM particle mass
    *fdm_tol = 1e-2;
    *verbose = 1;

    eosBM->type == TAB;
    strcpy(eosBM->filepath, "/path/to/eos/file");
    eosBM->n_P = 1;
    eosDM->type == TAB;
    strcpy(eosDM->filepath, "/path/to/eos/file");
    eosDM->n_P = 1;

    star_props->e_center = 2;

    // Delfault diff. rot params
    DiffRotBM->A = 0.8;
    DiffRotBM->B = 0.6;
    DiffRotBM->lambda1 = 2.0;
    DiffRotBM->lambda2 = 0.5;
    DiffRotBM->p = 1;
    DiffRotBM->csi = 3;

    DiffRotDM->A = 0.8;
    DiffRotDM->B = 0.6;
    DiffRotDM->lambda1 = 2.0;
    DiffRotDM->lambda2 = 0.5;
    DiffRotDM->p = 1;
    DiffRotDM->csi = 3;

    strcpy(input_filename, "config.d");
    strcpy(output2d_filename, "output.h5");
    strcpy(outputS_filename, "single.dat");
    strcpy(output_name, "auto");
    strcpy(id_file, "0");
}

void outfiles_name(char *output2d_filename, char *outputS_filename, char *output_name, size_t bufsize) {
    snprintf(outputS_filename, bufsize, "%s.dat", output_name);
    snprintf(output2d_filename, bufsize, "%s.h5", output_name);
}

void zero_DM(struct evolution_variables *ev, struct stellar_properties *DM_props, double *r_ratio_DM_target, double *r_e_DM) {
    DM_props->e_center = 1e-3;
    DM_props->RotType = UNIFORM;
    *r_ratio_DM_target = 1;

    *r_e_DM = 0.0;
    for (int s = 1; s <= SDIV; s++) {
        for (int m = 1; m <= MDIV; m++) {
            ev->energyDM[s][m] = 0.0;
            ev->pressureDM[s][m] = 0.0;
            ev->velocity_sqDM[s][m] = 0.0;
            ev->enthalpyDM[s][m] = 0.0;
            ev->Omega_diffDM[s][m] = 0.0;
        }
    }
}

void copy_fields(struct evolution_variables *Dest, struct evolution_variables *Src) {
    int s, m;

    for (s = 1; s <= SDIV; s++) {
        for (m = 1; m <= MDIV; m++) {
            Dest->rho[s][m] = Src->rho[s][m];
            Dest->gama[s][m] = Src->gama[s][m];
            Dest->alpha[s][m] = Src->alpha[s][m];
            Dest->omega[s][m] = Src->omega[s][m];
            Dest->Omega_diffBM[s][m] = Src->Omega_diffBM[s][m];
            Dest->energyBM[s][m] = Src->energyBM[s][m];
            Dest->pressureBM[s][m] = Src->pressureBM[s][m];
            Dest->enthalpyBM[s][m] = Src->enthalpyBM[s][m];
            Dest->velocity_sqBM[s][m] = Src->velocity_sqBM[s][m];
            Dest->Omega_diffDM[s][m] = Src->Omega_diffDM[s][m];
            Dest->energyDM[s][m] = Src->energyDM[s][m];
            Dest->pressureDM[s][m] = Src->pressureDM[s][m];
            Dest->enthalpyDM[s][m] = Src->enthalpyDM[s][m];
            Dest->velocity_sqDM[s][m] = Src->velocity_sqDM[s][m];
        }
    }
}

double bilinear_interpolator(double x[], double y[], double **f, int nx, int ny, double x0, double y0) {
    int ix = nx / 2, iy = ny / 2;

    hunt(x, nx, x0, &ix);
    hunt(y, ny, y0, &iy);
    ix += 1;
    iy += 1;

    double x1 = x[ix], x2 = x[ix + 1];
    double y1 = y[iy], y2 = y[iy + 1];

    double f11 = f[ix][iy];
    double f21 = f[ix + 1][iy];
    double f12 = f[ix][iy + 1];
    double f22 = f[ix + 1][iy + 1];

    double denom_x = x2 - x1;
    double denom_y = y2 - y1;
    double wx1 = (x2 - x0) / denom_x;
    double wx2 = (x0 - x1) / denom_x;
    double wy1 = (y2 - y0) / denom_y;
    double wy2 = (y0 - y1) / denom_y;

    return f11 * wx1 * wy1 +
           f21 * wx2 * wy1 +
           f12 * wx1 * wy2 +
           f22 * wx2 * wy2;
}

void regrid(double *s_gp, double *tmp_sgp, int tmp_sdiv, double *mu, double *tmp_mu, int tmp_mdiv, struct evolution_variables *ev, struct evolution_variables *tmp_ev) {
    if ((tmp_sdiv == SDIV) && (tmp_mdiv == MDIV)) {
        copy_fields(ev, tmp_ev);
        return;
    }

    DPRINT("Regridding: id is (%d, %d) while desired output is (%d, %d)\n", tmp_sdiv, tmp_mdiv, SDIV, MDIV);

    for (int s = 1; s <= SDIV; s++) {
        for (int m = 1; m <= MDIV; m++) {
            ev->rho[s][m] = bilinear_interpolator(tmp_sgp, tmp_mu, tmp_ev->rho, tmp_sdiv, tmp_mdiv, s_gp[s], mu[m]);
            ev->gama[s][m] = bilinear_interpolator(tmp_sgp, tmp_mu, tmp_ev->gama, tmp_sdiv, tmp_mdiv, s_gp[s], mu[m]);
            ev->alpha[s][m] = bilinear_interpolator(tmp_sgp, tmp_mu, tmp_ev->alpha, tmp_sdiv, tmp_mdiv, s_gp[s], mu[m]);
            ev->omega[s][m] = bilinear_interpolator(tmp_sgp, tmp_mu, tmp_ev->omega, tmp_sdiv, tmp_mdiv, s_gp[s], mu[m]);
            ev->Omega_diffBM[s][m] = bilinear_interpolator(tmp_sgp, tmp_mu, tmp_ev->Omega_diffBM, tmp_sdiv, tmp_mdiv, s_gp[s], mu[m]);
            ev->energyBM[s][m] = bilinear_interpolator(tmp_sgp, tmp_mu, tmp_ev->energyBM, tmp_sdiv, tmp_mdiv, s_gp[s], mu[m]);
            ev->pressureBM[s][m] = bilinear_interpolator(tmp_sgp, tmp_mu, tmp_ev->pressureBM, tmp_sdiv, tmp_mdiv, s_gp[s], mu[m]);
            ev->enthalpyBM[s][m] = bilinear_interpolator(tmp_sgp, tmp_mu, tmp_ev->enthalpyBM, tmp_sdiv, tmp_mdiv, s_gp[s], mu[m]);
            ev->velocity_sqBM[s][m] = bilinear_interpolator(tmp_sgp, tmp_mu, tmp_ev->velocity_sqBM, tmp_sdiv, tmp_mdiv, s_gp[s], mu[m]);
            ev->Omega_diffDM[s][m] = bilinear_interpolator(tmp_sgp, tmp_mu, tmp_ev->Omega_diffDM, tmp_sdiv, tmp_mdiv, s_gp[s], mu[m]);
            ev->energyDM[s][m] = bilinear_interpolator(tmp_sgp, tmp_mu, tmp_ev->energyDM, tmp_sdiv, tmp_mdiv, s_gp[s], mu[m]);
            ev->pressureDM[s][m] = bilinear_interpolator(tmp_sgp, tmp_mu, tmp_ev->pressureDM, tmp_sdiv, tmp_mdiv, s_gp[s], mu[m]);
            ev->enthalpyDM[s][m] = bilinear_interpolator(tmp_sgp, tmp_mu, tmp_ev->enthalpyDM, tmp_sdiv, tmp_mdiv, s_gp[s], mu[m]);
            ev->velocity_sqDM[s][m] = bilinear_interpolator(tmp_sgp, tmp_mu, tmp_ev->velocity_sqDM, tmp_sdiv, tmp_mdiv, s_gp[s], mu[m]);
        }
    }
}

void load_initial_data(char *id_file, double *s_gp, double *mu, struct evolution_variables *evolution_functions, struct evolution_variables *tmp_func, struct stellar_properties *star_props, struct stellar_properties *DM_props, struct EOS *eosBM, struct EOS *eosDM, double *r_ratio_BM, double *r_ratio_DM, double *Omega_BM, double *Omega_DM, double *r_e_BM, double *r_e_DM, double *s_e_BM, double *s_e_DM, struct DiffRotParams *DiffRotBM, struct DiffRotParams *DiffRotDM, bool counter, double *Mass_BM, double *Mass_0_BM, double *J_BM, double *R_e_BM, double *v_plus, double *v_minus, double *Omega_K_BM, double *Mass_DM, double *Mass_0_DM, double *J_DM, double *R_e_DM, double MDM, double *R_p_BM, double *R_p_DM, double *Omega_K_DM, double *T_BM, double *T_DM, double *W_BM, double *W_DM, double *GRV2, double *GRV3, int verbose) {
    int tmp_sdiv, tmp_mdiv;
    double tmp_sgp[16000], tmp_mu[16000];

    hdf5_read_var(tmp_sgp, &tmp_sdiv, tmp_mu, &tmp_mdiv, id_file, tmp_func, star_props, DM_props, eosBM, eosDM, r_ratio_BM, r_ratio_DM, Omega_BM, Omega_DM, r_e_BM, r_e_DM, DiffRotBM, DiffRotDM);

    regrid(s_gp, tmp_sgp, tmp_sdiv, mu, tmp_mu, tmp_mdiv, evolution_functions, tmp_func);

    make_center(eosBM, star_props->e_center, &(star_props->p_center), &(star_props->h_center));
    make_center(eosDM, DM_props->e_center, &(DM_props->p_center), &(DM_props->h_center));

    // mass_radius(s_gp, mu, eosBM, evolution_functions, *r_ratio_BM, *r_e_BM, *Omega_BM, Mass_BM, Mass_0_BM, J_BM, R_e_BM, v_plus, v_minus, Omega_K_BM, star_props, eosDM, *r_ratio_DM, *r_e_DM, *Omega_DM, Mass_DM, Mass_0_DM, J_DM, R_e_DM, DM_props, MDM, R_p_BM, R_p_DM, Omega_K_DM, T_BM, T_DM, W_BM, W_DM, counter, GRV2, GRV3);

    // double fdm = *Mass_DM / (*Mass_BM + *Mass_DM);
    // double mid_JDM = C * *J_DM / G / SQ(MSUN);
    // double Jtot = C * (*J_BM + *J_DM) / G / SQ(MSUN);

    // if (verbose >= STANDARD) {
    //     printf("Loaded data from %s\n", id_file);
    //     printf("Summary of loaded data:\n");
    //     printf("GRV2 = %g  (Eq. 61 of Astron. Astrophys. Suppl. Ser. 132, 431-454) \n", *GRV2);
    //     printf("GRV3 = %g  (Eq. 62 of Astron. Astrophys. Suppl. Ser. 132, 431-454) \n", *GRV3);
    //     printf("Dark Matter fraction: %g\n", fdm);
    //     printf("Total angular momentum: J_B + J_DM = %g + %g = %g\n\n", Jtot - mid_JDM, mid_JDM, Jtot);
    //     printf("Baryonic Matter:\n");
    //     printf("  %-15s %-15s %-15s %-15s %-15s %-15s\n", "Mass (Msun)", "Mass_0 (Msun)", "R_e (Km)", "J / M^2", "Omega BM (Hz)", "Omega_K BM (Hz)");
    //     printf("  %-15.8g %-15.8g %-15.8f %-15.8f %-15.8f %-15.8f\n",
    //            *Mass_BM / MSUN, *Mass_0_BM / MSUN, *R_e_BM / 1.0e5, (C * *J_BM / (G * *Mass_BM * *Mass_BM)), star_props->RotType == UNIFORM ? evolution_functions->Omega_diffBM[1][1] * C / sqrt(KAPPA) : NAN, *Omega_K_BM);
    //     printf("\nDark Matter:\n");
    //     printf("  %-15s %-15s %-15s %-15s %-15s %-15s\n", "Mass (Msun)", "Mass_0 (Msun)", "R_e (Km)", "J / M^2", "Omega DM (Hz)", "Omega_K DM (Hz)");
    //     printf("  %-15.8g %-15.8g %-15.8f %-15.8f %-15.8f %-15.8f\n",
    //            *Mass_DM / MSUN, *Mass_0_DM / MSUN, *R_e_DM / 1.0e5, (C * *J_DM / (G * *Mass_DM * *Mass_DM + DBL_EPSILON)), DM_props->RotType == UNIFORM ? evolution_functions->Omega_diffDM[1][1] * C / sqrt(KAPPA) : NAN, *Omega_K_DM);
    // }
    
    double r_e = fmax(*r_e_BM, *r_e_DM);
    *s_e_BM = *r_e_BM / (*r_e_BM + r_e);
    *s_e_DM = *r_e_DM / (*r_e_DM + r_e);
}

void do_jconstant(double *s_gp, double *mu, struct evolution_variables *ev, struct EOS *eosBM, struct EOS *eosDM, struct stellar_properties *star_props, struct stellar_properties *DM_props, struct DiffRotParams *DiffRotBM, struct DiffRotParams *DiffRotDM, double old_rrBM, double old_rrDM, double *r_e_BM, double *Omega_BM, double *s_e_BM, double *r_e_DM, double *Omega_DM, double *s_e_DM, bool *outOfiter, double *Omega_e_BM, double *Omega_e_DM, double r_ratio_BM, double r_ratio_DM, bool *BM_first_diff_done, bool *DM_first_diff_done, bool zero, double accuracy, double cf, bool counter, int verbose) {
    bool doSpecialSpin = false;
    int BMog_RotType, DMog_RotType;
    double A_BM_orig, A_DM_orig;
    int n_nearest = SDIV / 2;

    if ((star_props->RotType != UNIFORM && star_props->RotType != DIFF) &&
        (old_rrBM == 1 && r_ratio_BM != 1) &&
        !*BM_first_diff_done) {
        BMog_RotType = star_props->RotType;
        star_props->RotType = DIFF;
        A_BM_orig = DiffRotBM->A;
        DiffRotBM->A = 1;
        *BM_first_diff_done = true;
        doSpecialSpin = true;
    }

    DMog_RotType = DM_props->RotType;
    if ((DM_props->RotType != UNIFORM && DM_props->RotType != DIFF) &&
        (old_rrDM == 1 && r_ratio_DM != 1) &&
        !*DM_first_diff_done) {
        A_DM_orig = DiffRotDM->A;
        DiffRotDM->A = 1;
        DM_props->RotType = DIFF;
        *DM_first_diff_done = true;
        doSpecialSpin = true;
    }

    if (doSpecialSpin) {
        printf("Doing j-constant rotation for: %s %s\n", *BM_first_diff_done ? "BM" : "", *DM_first_diff_done ? "DM" : "");
        spin(s_gp, mu, eosBM, star_props, ev, accuracy, cf,
             r_ratio_BM, r_e_BM, Omega_BM, *s_e_BM,
             eosDM, DM_props, r_ratio_DM, r_e_DM, Omega_DM, *s_e_DM,
             verbose, outOfiter, counter,
             Omega_e_BM, DiffRotBM, Omega_e_DM, DiffRotDM, zero);

        printf("Applying parabola to Omega BM\n");
        int sLimit = r_e_BM > r_e_DM ? 300 : 60;

        double Omega_diff_mu_0[SDIV / 2];
        for (int s = 1; s <= SDIV / 2 - 1; s++) {
            Omega_diff_mu_0[s] = ev->Omega_diffBM[s][1];
        }
        double s_max = interp(Omega_diff_mu_0, s_gp, SDIV / 2 - 1, ev->Omega_diffBM[sLimit][1], &n_nearest);
        for (int s = 1; s <= sLimit; s++) {
            double s_tmp = interp(Omega_diff_mu_0, s_gp, SDIV / 2 - 1, ev->Omega_diffBM[s][1], &n_nearest);

            for (int m = 1; m <= MDIV; m++) {
                ev->Omega_diffBM[s][m] *= (-0.25 * SQ(s_tmp - s_max) / SQ(s_max) + 1);
            }
        }

        if (!zero && (DMog_RotType != UNIFORM && DMog_RotType != DIFF)) {
            printf("Applying parabola to Omega DM\n");
            sLimit = r_e_DM > r_e_BM ? 60 : 300;
            for (int s = 1; s <= SDIV / 2 - 1; s++) {
                Omega_diff_mu_0[s] = ev->Omega_diffDM[s][1];
            }
            s_max = interp(Omega_diff_mu_0, s_gp, SDIV / 2 - 1, ev->Omega_diffDM[sLimit][1], &n_nearest);
            for (int s = 1; s <= sLimit; s++) {
                double s_tmp = interp(Omega_diff_mu_0, s_gp, SDIV / 2 - 1, ev->Omega_diffDM[s][1], &n_nearest);

                for (int m = 1; m <= MDIV; m++) {
                    ev->Omega_diffDM[s][m] *= (-0.25 * SQ(s_tmp - s_max) / SQ(s_max) + 1);
                }
            }
        }

        if (verbose >= DEBUG) {
            FILE *fptro;
            fptro = fopen("omegaBM_parab.dat", "w");
            int m = 1;
            for (int s = 1; s <= 2 * SDIV / 3; s++) {
                fprintf(fptro, "%f %f %e\n", s_gp[s], mu[m], ev->Omega_diffBM[s][m]);
            }
            fclose(fptro);

            fptro = fopen("omegaDM_parab.dat", "w");
            for (int s = 1; s <= 2 * SDIV / 3; s++) {
                fprintf(fptro, "%f %f %e\n", s_gp[s], mu[m], ev->Omega_diffDM[s][m]);
            }
            fclose(fptro);
        }

        if (*BM_first_diff_done) {
            star_props->RotType = BMog_RotType;
            DiffRotBM->A = A_BM_orig;
        }
        if (*DM_first_diff_done) {
            DM_props->RotType = DMog_RotType;
            DiffRotDM->A = A_DM_orig;
        }
    }
}