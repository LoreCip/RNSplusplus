/**************************************************************************
 *
 *	          ___           ___           ___
 *	         /\  \         /\__\         /\  \
 *	        /::\  \       /::|  |       /::\  \
 *	       /:/\:\  \     /:|:|  |      /:/\ \  \
 *	      /::\~\:\  \   /:/|:|  |__   _\:\~\ \  \
 *	     /:/\:\ \:\__\ /:/ |:| /\__\ /\ \:\ \ \__\
 *	     \/_|::\/:/  / \/__|:|/:/  / \:\ \:\ \/__/
 *	        |:|::/  /      |:/:/  /   \:\ \:\__\
 *	        |:|\/__/       |::/  /     \:\/:/  /
 *	        |:|  |         /:/  /       \::/  /
 *	         \|__|         \/__/         \/__/
 *
 /**************************************************************************
 *
 * RNS_Diff_Jbm.c
 *
 * Author(s): Originally by N. Stergioulas, modified by L. Cipriani
 *
 * This program computes axisymmetric equilibrium models of rotating neutron
 * stars composed of baryonic matter and dark matter, under general relativity
 * with differential rotation. It is tailored to construct configurations
 * with a specified **baryonic angular momentum (J_BM)** and a **target dark
 * matter angular momentum fraction (fJ_DM = J_DM / J_tot)**.
 *
 * The code solves for equilibrium by varying the axis ratios and the dark
 * matter central energy density to satisfy constraints on the global dark
 * matter fraction and angular momenta. It adjusts `r_ratio_BM` and evolves
 * toward a configuration satisfying:
 *      - J_BM = user-specified value
 *      - J_DM = fJ_DM × (J_BM + J_DM)
 *
 * It supports loading previous configurations from HDF5 to improve convergence,
 * and outputs include bulk properties, virial diagnostics (GRV2/GRV3), and
 * optional 2D data dumps.
 *
 *
 **************************************************************************/
#define _GNU_SOURCE
#include <fenv.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "consts.h"
#include "equil.h"
#include "equil_util.h"
#include "main_util.h"
#include "nrutil.h"
#include "output.h"
#include "parfile.h"

int main(int argc, char **argv) {
    // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

    // EOS and stellar properties
    struct EOS eosBM;
    struct EOS eosDM;
    struct stellar_properties star_props, tmp_starprops;
    struct stellar_properties DM_props, tmp_DMprops;

    // Evolution variables
    struct evolution_variables evolution_functions, tmp_func;
    allocate_evolution_variables(&evolution_functions);
    allocate_evolution_variables(&tmp_func);

    // Differential rotation parameters
    struct DiffRotParams DiffRotBM, DiffRotDM;

    // Grid arrays
    double s_gp[SDIV + 1];  // s grid points
    double mu[MDIV + 1];    // mu grid points

    // Axis ratios and related targets/steps
    double r_ratio_BM = 1, r_ratio_DM = 1;
    double r_ratio_BM_target, r_ratio_DM_target;
    double r_ratio_step;

    // Equatorial and surface radii
    double r_e_BM, r_e_DM;
    double s_e_BM, s_e_DM;

    // Masses and angular momenta
    double Mass_BM, Mass_DM;
    double Mass_0_BM, Mass_0_DM;
    double J_BM, J_DM, JBM_Target, JDM_Target, Jtot_Target, fJ_DM;

    // Radii and rotation
    double R_e_BM, R_e_DM, R_p_BM, R_p_DM;
    double Omega_BM, Omega_DM;
    double Omega_K_BM, Omega_K_DM;
    double Omega_e_BM = 0, Omega_e_DM = 0;

    // Rotational velocities
    double *v_plus;   // velocity of co-rotating particle wrt ZAMO
    double *v_minus;  // velocity of counter-rotating particle wrt ZAMO
    v_plus = dvector(1, SDIV);
    v_minus = dvector(1, SDIV);

    // Dark matter fraction and related
    double fdm, fdm_target, fdm_tol, MDM;
    bool zero;

    // Energetics
    double T_BM, T_DM, W_BM, W_DM;

    // General relativistic virial identities
    double GRV2, GRV3;

    // Numerical parameters
    double cf, accuracy;
    bool load_id;

    // Output flags and verbosity
    bool OUTPUT2D, OUTPUT1D;
    int verbose;
    bool counter = false;
    bool outOfiter = false;  // Out-of-iteration flag

    // Filenames and output names
    size_t buffer_size = 512;
    char input_filename[buffer_size], output2d_filename[buffer_size], outputS_filename[buffer_size], output_name[buffer_size], id_file[buffer_size];

    /* SOME DEFAULT VALUES FOR WRITING THE PARAMETER FILE */
    set_default_values(&cf, &accuracy, &r_ratio_BM, &r_ratio_DM, &fdm_target, &MDM, &fdm_tol, &verbose, &star_props, &DM_props, &eosBM, &eosDM, &DiffRotBM, &DiffRotDM, input_filename, output2d_filename, outputS_filename, output_name, id_file);

    /*********************************/
    /*      PARAMETERS               */
    /*********************************/

    bool skip = false;
    int c;
    for (c = getopt(argc, argv, "+hc:j:d:");
         c != EOF;
         c = getopt(argc, argv, "+hc:j:d:")) {
        switch (c) {
            case 'c':
                sscanf(optarg, "%s", input_filename);
                skip = true;
                break;
            case 'j':
                sscanf(optarg, "%lf", &JBM_Target);
                break;
            case 'd':
                sscanf(optarg, "%lf", &fJ_DM);
                break;
            case 'h':
                printf("\n");
                printf("Help:\n");
                printf("Run one time to generate a template for the configuration file:\n");
                printf("\t>> ./RNS\n\n");
                printf("Run \n\t>> ./RNS -c config.d [-j JBM_Target] [-d fJ_DM]\n");
                printf("to start the computations.\n");
                return 0;
            default:
                fprintf(stderr, "Unknown option -%c\n", optopt);
                return 1;
        }
    }

    /* If no arguments, write config.d*/
    if (argc <= 2 && !skip) {
        writeConfig(input_filename, &eosBM, &star_props, accuracy, r_ratio_BM, OUTPUT1D, OUTPUT2D, &eosDM, &DM_props, r_ratio_DM, fdm_target, fdm_tol, MDM, verbose, counter, &DiffRotBM, &DiffRotDM, r_ratio_step);
        printf("Run \n\t>> ./RNS -c config.d\nto start the computations.\n");
        return 0;
    }
    readConfig(input_filename, &eosBM, &star_props, &accuracy, &r_ratio_BM_target, &OUTPUT1D, &OUTPUT2D, &eosDM, &DM_props, &r_ratio_DM_target, &fdm_target, &fdm_tol, &MDM, &verbose, &counter, &DiffRotBM, &DiffRotDM, &r_ratio_step, id_file, output_name, &load_id);

    outfiles_name(output2d_filename, outputS_filename, output_name, buffer_size);

    if (eosBM.type == TAB) {
        load_eos(&eosBM);
        star_props.e_surface = 7.8 * C * C * KSCALE;
        star_props.p_surface = 1.01e8 * KSCALE;
        star_props.enthalpy_min = 1.0 / (C * C);
    } else if (eosBM.type == POLY) {
        star_props.e_surface = 0.0;
        star_props.p_surface = 0.0;
        star_props.enthalpy_min = 0.0;
    }

    if (eosDM.type == TAB) {
        load_eos(&eosDM);
        DM_props.e_surface = pow(10., eosDM.log_e_tab[1]);
        DM_props.p_surface = pow(10., eosDM.log_p_tab[1]);
        DM_props.enthalpy_min = pow(10., eosDM.log_h_tab[1]);
    } else if (eosDM.type == POLY) {
        DM_props.e_surface = 0.0;
        DM_props.p_surface = 0.0;
        DM_props.enthalpy_min = 0.0;
    }
    double eDM_min = (eosDM.type == TAB) ? pow(10., eosDM.log_e_tab[1]) : 0.0;

    make_grid(s_gp, mu);

    make_center(&eosBM, star_props.e_center, &(star_props.p_center), &(star_props.h_center));
    make_center(&eosDM, DM_props.e_center, &(DM_props.p_center), &(DM_props.h_center));

    /* Here it computes a non rotating configuration as a first guess for the rns */
    double massBM, massDM;
    sphere(s_gp, &eosBM, &star_props, &evolution_functions, &r_e_BM, &s_e_BM, &eosDM, &DM_props, &r_e_DM, &s_e_DM, &massBM, &massDM);

    zero = false;
    if (fdm_target == 0) {
        zero = true;
        zero_DM(&evolution_functions, &DM_props, &r_ratio_DM_target, &r_e_DM);
    }

    if (load_id) {
        load_initial_data(id_file, s_gp, mu, &evolution_functions, &tmp_starprops, &DM_props, &eosBM, &eosDM, &r_ratio_BM, &r_ratio_DM, &Omega_BM, &Omega_DM, &r_e_BM, &r_e_DM, &s_e_BM, &s_e_DM, &DiffRotBM, &DiffRotDM, counter, &Mass_BM, &Mass_0_BM, &J_BM, &R_e_BM, v_plus, v_minus, &Omega_K_BM, &Mass_DM, &Mass_0_DM, &J_DM, &R_e_DM, MDM, &R_p_BM, &R_p_DM, &Omega_K_DM, &T_BM, &T_DM, &W_BM, &W_DM, &GRV2, &GRV3, verbose);

        // Bring ec_BM to desired value
        while (fabs(tmp_starprops.e_center - star_props.e_center) > 1e-8) {
            double diff = star_props.e_center - tmp_starprops.e_center;
            double step = 0.1;
            if (fabs(diff) < step) {
                tmp_starprops.e_center = star_props.e_center;
            } else {
                tmp_starprops.e_center += SIGN(1, diff) * step;
            }
            make_center(&eosBM, tmp_starprops.e_center, &(tmp_starprops.p_center), &(tmp_starprops.h_center));
            DPRINT("Computing ec_BM = %g\n", tmp_starprops.e_center);

            spin(s_gp, mu, &eosBM, &tmp_starprops, &evolution_functions, accuracy, cf, r_ratio_BM, &r_e_BM, &Omega_BM, s_e_BM, &eosDM, &DM_props, r_ratio_DM, &r_e_DM, &Omega_DM, s_e_DM, verbose, &outOfiter, counter, &Omega_e_BM, &DiffRotBM, &Omega_e_DM, &DiffRotDM, zero);

            if (outOfiter) {
                printf("Too many iterations trying to reche desired BM e_center. Exiting.\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    if (verbose >= STANDARD) {
        if (eosBM.type == TAB)
            printf("Chosen baryonic EOS file: %s\n", eosBM.filepath);
        else if (eosBM.type == POLY)
            printf("Baryonic polytropic EoS with N=%f\n", eosBM.n_P);
        if (eosDM.type == TAB)
            printf("Chosen dark matter EOS file: %s\n", eosDM.filepath);
        else if (eosDM.type == POLY)
            printf("Dark matter polytropic EoS with N=%f\n", eosDM.n_P);
        printf("Grid size: MDIVxSDIV=%dx%d\n\n", MDIV, SDIV);

        printf("**********************************************************\n");
        printf("**********************************************************\n");
        printf("Computing model with central baryonic energy density: %g\n", star_props.e_center);
        printf("Target dark matter fraction: %g %%\n", fdm_target * 100);
        printf("Target axis ratio for Baryonic Matter: %.4g\n", r_ratio_BM_target);
        if (!zero) printf("Target axis ratio for Dark Matter    : %.4g\n", r_ratio_DM_target);
        printf("**********************************************************\n");
        printf("**********************************************************\n");
    }

    ///////////////////////////
    //  END OF COMMON LOGIC  //
    ///////////////////////////

    // --- Temporary variables for iterative solution ---
    double tmp_reBM, tmp_OmegaBM, tmp_seBM;               // Temporary BM radius, Omega, s_e
    double tmp_reDM, tmp_OmegaDM, tmp_seDM;               // Temporary DM radius, Omega, s_e
    double tmp_OmegaeBM, tmp_OmegaeDM;                    // Temporary equatorial Omegas
    double tmp_rratioDM;                                  // Temporary DM axis ratio
    double mid_JDM, mid_JDM_old;                          // Temporary DM angular momentum (and previous)
    bool tmp_BM_first_diff_done, tmp_DM_first_diff_done;  // Flags for first differential rotation step

    // --- Flags for first differential rotation step ---
    bool BM_first_diff_done = false;
    bool DM_first_diff_done = false;

    // --- DM central energy density and fraction tracking ---
    double DM_e_center_old = DM_props.e_center;
    double fdm_old;

    // --- Previous axis ratios for j-constant step ---
    double old_rrBM = 1.1, old_rrDM = 1.1;

    // --- Angular momentum bookkeeping ---
    double Jtot = 0;

    // --- Bracketing and step variables for BM axis ratio iteration ---
    double r_ratio_BM_hi = 1;         // Upper bound for BM axis ratio
    double J_rhi = 0;                 // J at upper bound
    double r_ratio_BM_lo = 0.0;       // Lower bound for BM axis ratio
    double J_rlo = 10000;             // J at lower bound
    bool accepted_last_step = false;  // Flag to detect oscillation

    // --- Step sizes for axis ratio iteration ---
    double r_ratio_step_BM = r_ratio_step;
    double r_ratio_step_DM = r_ratio_step;

    // --- Timing variables ---
    struct timespec start, end;

    // TARGETS //
    Jtot_Target = JBM_Target / (1 - fJ_DM);
    JDM_Target = fJ_DM * Jtot_Target;

    // Task: find ec_DM and r_ratio_BM that give fdm and JBM for a fixed ec_BM and r_ratio_DM
    // Problem: r_ratio_BM must be always decreasing and if r_ratio_dm is small BM will not converge

    int Jiter = 0;
    // Loop over r_ratio_BM
    while (Jiter < MAXIT) {
        DPRINT("J_BM loop, iteration = %d\n", Jiter);

        if (verbose >= STANDARD) {
            printf("==========================================================\n");
            printf("\nIteration %d...\n", Jiter);
            printf("r_ratio_BM = %g      r_ratio_DM = %.4g\n", r_ratio_BM, r_ratio_DM);
            printf("==========================================================\n");
        }

        DPRINT("Storing temporary values for BM and DM before evolution\n");

        // Copy the good previous solution in a tmp array and make it evolve
        copy_fields(&tmp_func, &evolution_functions);
        tmp_reBM = r_e_BM;
        tmp_OmegaBM = Omega_BM;
        tmp_seBM = s_e_BM;
        tmp_reDM = r_e_DM;
        tmp_OmegaDM = Omega_DM;
        tmp_seDM = s_e_DM;
        tmp_OmegaeBM = Omega_e_BM;
        tmp_OmegaeDM = Omega_e_DM;
        tmp_rratioDM = r_ratio_DM;
        tmp_BM_first_diff_done = BM_first_diff_done;
        tmp_DM_first_diff_done = DM_first_diff_done;

        tmp_DMprops.e_center = DM_props.e_center;
        tmp_DMprops.p_center = DM_props.p_center;
        tmp_DMprops.h_center = DM_props.h_center;
        tmp_DMprops.e_surface = DM_props.e_surface;
        tmp_DMprops.p_surface = DM_props.p_surface;
        tmp_DMprops.enthalpy_min = DM_props.enthalpy_min;
        tmp_DMprops.RotType = DM_props.RotType;

        // Loop to change r_ratio_DM and keep J_DM constant
        int Jdm_iter = 0;
        r_ratio_step_DM = r_ratio_step;

        // Loop over DM central energy
        double r_ratio_DM_lo = 0.0, r_ratio_DM_hi = 1.0;
        double J_rlo_DM = 10000, J_rhi_DM = 0.0;
        bool accepted_dm_last_step = false;
        while (Jdm_iter < MAXIT) {
            DPRINT("J_DM loop, iteration = %d\n", Jdm_iter);
            int iter = 0;
            while (iter < MAXIT) {
                DPRINT("f_DM loop, iteration = %d\n", iter);
                if (verbose >= STANDARD) {
                    printf("----------------------------------------------------------\n");
                    printf(" ec_BM = %g      ec_DM = %g\n", star_props.e_center, tmp_DMprops.e_center);
                    printf("----------------------------------------------------------\n");
                }
                outOfiter = false;

                // Check if J-constant first step is necessary
                // do_jconstant(s_gp, mu, &tmp_func, &eosBM, &eosDM, &star_props, &tmp_DMprops, &DiffRotBM, &DiffRotDM, old_rrBM, old_rrDM, &tmp_reBM, &tmp_OmegaBM, &tmp_seBM, &tmp_reDM, &tmp_OmegaDM, &tmp_seDM, &outOfiter, &tmp_OmegaeBM, &tmp_OmegaeDM, r_ratio_BM, tmp_rratioDM, &tmp_BM_first_diff_done, &tmp_DM_first_diff_done, zero, accuracy, cf, counter, verbose);

                // Evolve the temp functions
                DPRINT("Calling spin() with r_ratio_BM = %g, r_ratio_DM = %g\n", r_ratio_BM, tmp_rratioDM);
                clock_gettime(CLOCK_MONOTONIC, &start);
                spin(s_gp, mu, &eosBM, &star_props, &tmp_func, accuracy, cf, r_ratio_BM, &tmp_reBM, &tmp_OmegaBM, tmp_seBM, &eosDM, &tmp_DMprops, tmp_rratioDM, &tmp_reDM, &tmp_OmegaDM, tmp_seDM, verbose, &outOfiter, counter, &tmp_OmegaeBM, &DiffRotBM, &tmp_OmegaeDM, &DiffRotDM, zero);
                clock_gettime(CLOCK_MONOTONIC, &end);

                // Compute integral quantities
                mass_radius(s_gp, mu, &eosBM, &tmp_func, r_ratio_BM, tmp_reBM, tmp_OmegaBM, &Mass_BM, &Mass_0_BM, &J_BM, &R_e_BM, v_plus, v_minus, &Omega_K_BM, &star_props, &eosDM, tmp_rratioDM, tmp_reDM, tmp_OmegaDM, &Mass_DM, &Mass_0_DM, &J_DM, &R_e_DM, &tmp_DMprops, MDM, &R_p_BM, &R_p_DM, &Omega_K_DM, &T_BM, &T_DM, &W_BM, &W_DM, counter, &GRV2, &GRV3);

                fdm_old = fdm;
                fdm = Mass_DM / (Mass_BM + Mass_DM);
                mid_JDM_old = mid_JDM;
                mid_JDM = C * J_DM / G / SQ(MSUN);
                Jtot = C * (J_BM + J_DM) / G / SQ(MSUN);

                if (verbose >= STANDARD) {
                    double time_spent = (end.tv_sec - start.tv_sec) +
                                        (end.tv_nsec - start.tv_nsec) / 1e9;
                    printf("Time taken: %g seconds\n", time_spent);
                    printf("\nNew BM: A = %g   B = %g\n", DiffRotBM.A, DiffRotBM.B);
                    printf("New DM: A = %g   B = %g\n", DiffRotDM.A, DiffRotDM.B);
                    printf("\nGRV2 = %g  (Eq. 61 of Astron. Astrophys. Suppl. Ser. 132, 431-454) \n", GRV2);
                    printf("GRV3 = %g  (Eq. 62 of Astron. Astrophys. Suppl. Ser. 132, 431-454) \n", GRV3);
                    printf("\nFound Dark Matter fraction: %g\n", fdm);
                    printf("Found total angular momentum: J_B + J_DM = %g + %g = %g\n\n", Jtot - mid_JDM, mid_JDM, Jtot);
                }

                DPRINT("Computed fdm = %.6f (old = %.6f), mid_JDM = %.6f, Jtot = %.6f\n", fdm, fdm_old, mid_JDM, Jtot);

                if ((zero && fdm <= fdm_target) || (fabs(fdm - fdm_target) / fdm_target < fdm_tol)) {
                    DPRINT("DM fdm converged: fdm = %.6f, fdm_target = %.6f\n", fdm, fdm_target);
                    break;
                }

                double tmp = tmp_DMprops.e_center;
                DPRINT("Before update: e_center = %.6e, fdm = %.6f, fdm_old = %.6f\n", tmp, fdm, fdm_old);

                if (iter == 0) {
                    double scale = fabs(fdm - fdm_target) < 1e-1 ? 0.05 : 0.25;
                    double factor = 1 + SIGN(1, fdm_target - fdm) * scale;
                    double proposed = tmp_DMprops.e_center * factor;
                    if (proposed < eDM_min) {
                        // Bisection between min and current
                        tmp_DMprops.e_center = 0.5 * (tmp_DMprops.e_center + eDM_min);
                        DPRINT("  Initial step: bisection to min allowed e_center → %.6e\n", tmp_DMprops.e_center);
                    } else {
                        tmp_DMprops.e_center = proposed;
                        DPRINT("  Initial step: multiplying e_center by %.3f → %.6e\n", factor, tmp_DMprops.e_center);
                    }
                } else {
                    double delta = -(fdm - fdm_target) * (tmp_DMprops.e_center - DM_e_center_old) / (fdm - fdm_old);
                    double proposed = tmp_DMprops.e_center + delta;
                    if (proposed < eDM_min) {
                        // Bisection between min and current
                        tmp_DMprops.e_center = 0.5 * (tmp_DMprops.e_center + eDM_min);
                        DPRINT("  Secant step below min: bisection to min allowed e_center → %.6e\n", tmp_DMprops.e_center);
                    } else {
                        tmp_DMprops.e_center = proposed;
                        DPRINT("  Secant step: new e_center = %.6e\n", tmp_DMprops.e_center);
                    }
                }
                DM_e_center_old = tmp;

                iter++;
                make_center(&eosDM, tmp_DMprops.e_center, &(tmp_DMprops.p_center), &(tmp_DMprops.h_center));
            }
            DPRINT("Post f_DM inner loop: mid_JDM = %.6f, JDM_Target = %.6f\n", mid_JDM, JDM_Target);

            if (zero || fabs(mid_JDM - JDM_Target) / JDM_Target <= fdm_tol) {
                DPRINT("J_DM converged.\n");
                break;
            }

            // Check if JDM is below JDM_target
            if (mid_JDM <= JDM_Target) {
                // Accept and update upper bound
                r_ratio_DM_hi = tmp_rratioDM;
                J_rhi_DM = mid_JDM;

                DPRINT("r_ratio_DM_lo = %g J_rlo_DM = %g  r_ratio_DM_hi = %g  J_rhi_DM = %g\n", r_ratio_DM_lo, J_rlo_DM, r_ratio_DM_hi, J_rhi_DM);
                if (r_ratio_DM_lo != 0.0 && fabs(J_rhi_DM - J_rlo_DM) > 1e-12) {
                    double new_rr = r_ratio_DM_hi - (J_rhi_DM - JDM_Target) * (r_ratio_DM_hi - r_ratio_DM_lo) / (J_rhi_DM - J_rlo_DM);
                    tmp_rratioDM = fmax(fmin(new_rr, r_ratio_DM_hi), r_ratio_DM_lo);
                    DPRINT("DM accepted: secant refinement → r_ratio_DM = %g\n", tmp_rratioDM);
                } else {
                    tmp_rratioDM -= r_ratio_step_DM;  // Jdm_iter == 0 ? 0.005 : r_ratio_step_DM;
                    DPRINT("DM accepted: fallback decrement → r_ratio_DM = %g\n", tmp_rratioDM);
                }

                accepted_dm_last_step = true;
            } else {
                // Reject and update lower bound
                r_ratio_DM_lo = tmp_rratioDM;
                J_rlo_DM = mid_JDM;

                if (accepted_dm_last_step) {
                    r_ratio_step_DM /= 2;
                    accepted_dm_last_step = false;
                }

                DPRINT("r_ratio_DM_lo = %g J_rlo_DM = %g  r_ratio_DM_hi = %g  J_rhi_DM = %g\n", r_ratio_DM_lo, J_rlo_DM, r_ratio_DM_hi, J_rhi_DM);
                if (r_ratio_DM_hi != 1 && fabs(J_rlo_DM - J_rhi_DM) > 1e-12) {
                    double new_rr = r_ratio_DM_hi - (J_rhi_DM - JDM_Target) * (r_ratio_DM_hi - r_ratio_DM_lo) / (J_rhi_DM - J_rlo_DM);
                    tmp_rratioDM = fmax(fmin(new_rr, r_ratio_DM_hi), r_ratio_DM_lo);
                    DPRINT("DM rejected: secant step → r_ratio_DM = %g\n", tmp_rratioDM);
                } else {
                    tmp_rratioDM += r_ratio_step_DM;
                    DPRINT("DM rejected: fallbakc increase → r_ratio_DM = %g\n", tmp_rratioDM);
                }
            }

            Jdm_iter++;
        }

        DPRINT("Checking total J convergence: Jtot = %.6f, Jtot_Target = %.6f\n", Jtot, Jtot_Target);
        // When J_DM and fdm target are reached,
        // Check convergence
        if (fabs(Jtot - Jtot_Target) / Jtot_Target <= fdm_tol) {
            DPRINT("Jtot converged. Copying accepted values and exiting loop.\n");
            // Step accepted, copy values back
            copy_fields(&evolution_functions, &tmp_func);
            r_e_BM = tmp_reBM;
            Omega_BM = tmp_OmegaBM;
            s_e_BM = tmp_seBM;
            r_e_DM = tmp_reDM;
            Omega_DM = tmp_OmegaDM;
            s_e_DM = tmp_seDM;
            Omega_e_BM = tmp_OmegaeBM;
            Omega_e_DM = tmp_OmegaeDM;
            r_ratio_DM = tmp_rratioDM;
            BM_first_diff_done = tmp_BM_first_diff_done;
            DM_first_diff_done = tmp_DM_first_diff_done;

            DM_props.e_center = tmp_DMprops.e_center;
            DM_props.p_center = tmp_DMprops.p_center;
            DM_props.h_center = tmp_DMprops.h_center;

            // And exit
            break;
        }

        // Check if Jtot is below Jtarget
        if (Jtot <= Jtot_Target) {
            // Step accepted, copy values back
            copy_fields(&evolution_functions, &tmp_func);
            r_e_BM = tmp_reBM;
            Omega_BM = tmp_OmegaBM;
            s_e_BM = tmp_seBM;
            r_e_DM = tmp_reDM;
            Omega_DM = tmp_OmegaDM;
            s_e_DM = tmp_seDM;
            Omega_e_BM = tmp_OmegaeBM;
            Omega_e_DM = tmp_OmegaeDM;
            r_ratio_DM = tmp_rratioDM;
            BM_first_diff_done = tmp_BM_first_diff_done;
            DM_first_diff_done = tmp_DM_first_diff_done;

            DM_props.e_center = tmp_DMprops.e_center;
            DM_props.p_center = tmp_DMprops.p_center;
            DM_props.h_center = tmp_DMprops.h_center;

            old_rrBM = r_ratio_BM;

            // Update upper bound and function value
            r_ratio_BM_hi = r_ratio_BM;
            J_rhi = Jtot;

            DPRINT("r_ratio_BM_lo = %g J_rlo_BM = %g  r_ratio_BM_hi = %g  J_rhi_BM = %g\n", r_ratio_BM_lo, J_rlo, r_ratio_BM_hi, J_rhi);
            double denom = (J_rhi - J_rlo);
            // To move r_ratio_BM, use fixed decrement until an estimate of J_rlo can be given
            // Once it is bracketed, try to use secant
            if (r_ratio_BM_lo != 0.0 && fabs(denom) > 1e-12) {
                double r_ratio_new = r_ratio_BM_hi - (J_rhi - Jtot_Target) * (r_ratio_BM_hi - r_ratio_BM_lo) / denom;
                // Clamp within bounds
                r_ratio_BM = fmax(fmin(r_ratio_new, r_ratio_BM_hi), r_ratio_BM_lo);
                DPRINT("BM accepted: secant refinement → r_ratio_BM = %g\n", r_ratio_BM);
            } else {
                r_ratio_BM -= r_ratio_step_BM;  // Jiter == 0 ? 0.005 : r_ratio_step_BM;
                DPRINT("BM accepted: fallback decrement → r_ratio_BM = %g\n", r_ratio_BM);
            }

            accepted_last_step = true;

            Jiter++;
            DPRINT("Step accepted, going to r_ratio_BM = %g\n", r_ratio_BM);

        } else {
            // If step is rejected, update lower bound for r_ratio_BM
            r_ratio_BM_lo = r_ratio_BM;
            J_rlo = Jtot;

            // Halven the step
            if (accepted_last_step) {
                r_ratio_step_BM /= 2;
                accepted_last_step = false;
            }

            DPRINT("r_ratio_BM_lo = %g J_rlo_BM = %g  r_ratio_BM_hi = %g  J_rhi_BM = %g\n", r_ratio_BM_lo, J_rlo, r_ratio_BM_hi, J_rhi);
            // Try to use secant method
            // If we are here we know for sure that the r_ratio_BM is between
            // r_ratio_BM_lo and r_ratio_BM_hi
            double denom = (J_rlo - J_rhi);
            if (r_ratio_BM_hi != 1.0 && fabs(denom) > 1e-12) {
                r_ratio_BM = r_ratio_BM_hi - (J_rhi - Jtot_Target) * (r_ratio_BM_hi - r_ratio_BM_lo) / denom;
                r_ratio_BM = fmax(fmin(r_ratio_BM, r_ratio_BM_hi), r_ratio_BM_lo);
                DPRINT("BM rejected: secant step → r_ratio_BM = %g\n", r_ratio_BM);
            } else {
                r_ratio_BM += r_ratio_step_BM;  // Jiter == 0 ? 0.005 : r_ratio_step_BM;
                DPRINT("BM accepted: fallback increment → r_ratio_BM = %g\n", r_ratio_BM);
            }
        }
        DPRINT("\n");
    }

    Omega_BM = star_props.RotType == UNIFORM ? evolution_functions.Omega_diffBM[1][1] * C / sqrt(KAPPA) : NAN;
    Omega_DM = DM_props.RotType == UNIFORM ? evolution_functions.Omega_diffDM[1][1] * C / sqrt(KAPPA) : NAN;
    if (verbose >= STANDARD) {
        printf("Baryonic Matter:\n");
        printf("  %-15s %-15s %-15s %-15s %-15s %-15s\n", "Mass (Msun)", "Mass_0 (Msun)", "R_e (Km)", "J / M^2", "Omega BM (Hz)", "Omega_K BM (Hz)");
        printf("  %-15.8g %-15.8g %-15.8f %-15.8f %-15.8f %-15.8f\n",
               Mass_BM / MSUN, Mass_0_BM / MSUN, R_e_BM / 1.0e5, (C * J_BM / (G * Mass_BM * Mass_BM)), Omega_BM, Omega_K_BM);
        if (!zero) {
            printf("\nDark Matter:\n");
            printf("  %-15s %-15s %-15s %-15s %-15s %-15s\n", "Mass (Msun)", "Mass_0 (Msun)", "R_e (Km)", "J / M^2", "Omega DM (Hz)", "Omega_K DM (Hz)");
            printf("  %-15.8g %-15.8g %-15.8f %-15.8f %-15.8f %-15.8f\n",
                   Mass_DM / MSUN, Mass_0_DM / MSUN, R_e_DM / 1.0e5, (C * J_DM / (G * Mass_DM * Mass_DM + DBL_EPSILON)), Omega_DM, Omega_K_DM);
        }
        printf("\n");
    }

    FILE *file;
    open_output(outputS_filename, &file);
    print_output(file, &star_props, R_e_BM, R_p_BM, Mass_BM, Mass_0_BM, &DM_props, R_e_DM, R_p_DM, Mass_DM, Mass_0_DM, fdm, Omega_BM, Omega_DM, J_BM, J_DM, r_ratio_BM, r_ratio_DM, Omega_K_BM, Omega_K_DM, T_BM, T_DM, W_BM, W_DM);
    close_output(file);

    /********************************************
        OUTPUTS
    *********************************************/

    if (OUTPUT2D) {
        char folder_path[256] = "./output";
        char folder_path_2d[256] = "./output/output_2d";
        char filename[512];

        // Create folders if they don't exist
        if (access(folder_path, F_OK) == -1) {
            if (mkdir(folder_path, 0777) == -1) {
                perror("Error creating folder");
                return 1;
            }
            printf("Folder created: %s\n", folder_path);
        }
        if (access(folder_path_2d, F_OK) == -1) {
            if (mkdir(folder_path_2d, 0777) == -1) {
                perror("Error creating folder for output_2d");
                return 1;
            }
            printf("Folder created: %s\n", folder_path_2d);
        }

        // Full path to output HDF5
        char output_path[1024];
        snprintf(output_path, sizeof(output_path), "%s/%s", folder_path_2d, output2d_filename);
        printf("Saving H5 file in: \n");
        printf("%s\n", output_path);

        // Save HDF5
        hdf5_save_var(s_gp, SDIV, mu, MDIV, output_path, &evolution_functions, &star_props,
                      &DM_props, &eosBM, &eosDM, &r_ratio_BM, &r_ratio_DM,
                      &Omega_BM, &Omega_DM, &r_e_BM, &r_e_DM,
                      &DiffRotBM, &DiffRotDM, input_filename);
    }

    return 0;
}