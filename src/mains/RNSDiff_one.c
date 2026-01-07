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
 * RNS_Diff_one.c
 *
 * Author(s): Originally by N. Stergioulas, modified by L. Cipriani
 *
 * This program computes axisymmetric, stationary equilibrium models of
 * rotating neutron stars with both baryonic and dark matter components
 * under general relativity, using the Komatsu-Eriguchi-Hachisu (KEH) method.
 *
 * It constructs equilibrium sequences for specified baryonic and dark
 * matter equations of state (EOS), solving for rotating configurations
 * with differential rotation profiles. Starting from a non-rotating
 * configuration, it iteratively decreases axis ratios to reach target
 * flattening and specified dark matter fraction.
 *
 * Output includes global properties (mass, radius, angular momentum,
 * rotational kinetic and gravitational binding energies), virial checks,
 * and optional 2D HDF5 snapshots of the metric and matter fields.
 *
 * Input parameters (EOS, central energy densities, rotation parameters, etc.)
 * are specified in a configuration file provided via command-line argument.
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

int main(int argc, char** argv) {
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
    double r_ratio_step = 0.05;

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
    double* v_plus;   // velocity of co-rotating particle wrt ZAMO
    double* v_minus;  // velocity of counter-rotating particle wrt ZAMO
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
    for (c = getopt(argc, argv, "+hc:");
         c != EOF;
         c = getopt(argc, argv, "+hc:")) {
        switch (c) {
            case 'c':
                sscanf(optarg, "%s", input_filename);
                skip = true;
                break;
            case 'h':
                printf("\n");
                printf("Help:\n");
                printf("Run one time to generate a template for the configuration file:\n");
                printf("\t>> ./RNS\n\n");
                printf("Run \n\t>> ./RNS -c config.d\nto start the computations.\n");
                return 0;
            default:
                break;
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
        struct DiffRotParams tmp_DiffRotBM, tmp_DiffRotDM;
        load_initial_data(id_file, s_gp, mu, &evolution_functions, &tmp_func, &tmp_starprops, &DM_props, &eosBM, &eosDM, &r_ratio_BM, &r_ratio_DM, &Omega_BM, &Omega_DM, &r_e_BM, &r_e_DM, &s_e_BM, &s_e_DM, &tmp_DiffRotBM, &tmp_DiffRotDM, counter, &Mass_BM, &Mass_0_BM, &J_BM, &R_e_BM, v_plus, v_minus, &Omega_K_BM, &Mass_DM, &Mass_0_DM, &J_DM, &R_e_DM, MDM, &R_p_BM, &R_p_DM, &Omega_K_DM, &T_BM, &T_DM, &W_BM, &W_DM, &GRV2, &GRV3, verbose);

        DiffRotBM.A = tmp_DiffRotBM.A;
        DiffRotBM.B = tmp_DiffRotBM.B;
        DiffRotDM.A = tmp_DiffRotDM.A;
        DiffRotDM.B = tmp_DiffRotDM.B;

        // Bring ec_BM to desired value
        while (fabs(tmp_starprops.e_center - star_props.e_center) > 1e-8) {
            double delta = star_props.e_center - tmp_starprops.e_center;
            if (fabs(delta) < 0.1) {
                tmp_starprops.e_center = star_props.e_center;
            } else {
                tmp_starprops.e_center += 0.1 * SIGN(1, delta);
            }
            make_center(&eosBM, tmp_starprops.e_center, &(tmp_starprops.p_center), &(tmp_starprops.h_center));
            DPRINT("Computing ec_BM = %g\n", tmp_starprops.e_center);

            spin(s_gp, mu, &eosBM, &tmp_starprops, &evolution_functions, accuracy, cf, r_ratio_BM, &r_e_BM, &Omega_BM, s_e_BM, &eosDM, &tmp_DMprops, r_ratio_DM, &r_e_DM, &Omega_DM, s_e_DM, verbose, &outOfiter, counter, &Omega_e_BM, &DiffRotBM, &Omega_e_DM, &DiffRotDM, zero);
            if (outOfiter) {
                printf("Too many iterations trying to reche desired BM e_center. Exiting.\n");
                exit(EXIT_FAILURE);
            }
        }

        // Bring ec_DM to desired value
        while (!zero && fabs(DM_props.e_center - tmp_DMprops.e_center) > 1e-8) {
            double delta = DM_props.e_center - tmp_DMprops.e_center;
            if (fabs(delta) < 0.1) {
                tmp_DMprops.e_center = DM_props.e_center;
            } else {
                tmp_DMprops.e_center += 0.1 * SIGN(1, delta);
            }

            make_center(&eosDM, tmp_DMprops.e_center, &(tmp_DMprops.p_center), &(tmp_DMprops.h_center));
            DPRINT("Computing ec_DM = %g\n", tmp_DMprops.e_center);

            spin(s_gp, mu, &eosBM, &tmp_starprops, &evolution_functions, accuracy, cf, r_ratio_BM, &r_e_BM, &Omega_BM, s_e_BM, &eosDM, &tmp_DMprops, r_ratio_DM, &r_e_DM, &Omega_DM, s_e_DM, verbose, &outOfiter, counter, &Omega_e_BM, &DiffRotBM, &Omega_e_DM, &DiffRotDM, zero);
            if (outOfiter) {
                printf("Too many iterations trying to reche desired DM e_center. Exiting.\n");
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
        printf("Target axis ratio for Dark Matter    : %.4g\n", r_ratio_DM_target);
        if (zero) printf("Ignoring Dark Matter component...\n");
        printf("**********************************************************\n");
        printf("**********************************************************\n");
    }

    ///////////////////////////
    //  END OF COMMON LOGIC  //
    ///////////////////////////

    // --- DM central energy density and fraction tracking ---
    double DM_e_center_old = DM_props.e_center;
    double fdm_old;

    // --- Timing variables ---
    struct timespec start, end;

    // Reach r_ratio target
    while (fabs(r_ratio_BM - r_ratio_BM_target) >= 1e-8 || fabs(r_ratio_DM - r_ratio_DM_target) >= 1e-8) {
        if (verbose >= STANDARD) {
            printf("----------------------------------------------------------\n");
            printf("r_ratio_BM = %.4g     r_ratio_DM = %.4g\n", r_ratio_BM, r_ratio_DM);
            printf("----------------------------------------------------------\n");
        }

        bool outOfiter = false;

        clock_gettime(CLOCK_MONOTONIC, &start);
        spin(s_gp, mu, &eosBM, &star_props, &evolution_functions, accuracy, cf, r_ratio_BM, &r_e_BM, &Omega_BM, s_e_BM, &eosDM, &DM_props, r_ratio_DM, &r_e_DM, &Omega_DM, s_e_DM, verbose, &outOfiter, counter, &Omega_e_BM, &DiffRotBM, &Omega_e_DM, &DiffRotDM, zero);
        clock_gettime(CLOCK_MONOTONIC, &end);

        if (outOfiter) {
            printf("This should not happen...");
        }

        if (fabs(r_ratio_BM - r_ratio_BM_target) < 1e-4 && fabs(r_ratio_DM - r_ratio_DM_target) < 1e-4) {
            break;
        }

        double delta = r_ratio_BM_target - r_ratio_BM;
        if (fabs(delta) < r_ratio_step) {
            r_ratio_BM = r_ratio_BM_target;
        } else {
            r_ratio_BM += r_ratio_step * SIGN(1, delta);
        }

        delta = r_ratio_DM_target - r_ratio_DM;
        if (fabs(delta) < r_ratio_step) {
            r_ratio_DM = r_ratio_DM_target;
        } else {
            r_ratio_DM += r_ratio_step * SIGN(1, delta);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &start);
    spin(s_gp, mu, &eosBM, &star_props, &evolution_functions, accuracy, cf, r_ratio_BM, &r_e_BM, &Omega_BM, s_e_BM, &eosDM, &DM_props, r_ratio_DM, &r_e_DM, &Omega_DM, s_e_DM, verbose, &outOfiter, counter, &Omega_e_BM, &DiffRotBM, &Omega_e_DM, &DiffRotDM, zero);
    clock_gettime(CLOCK_MONOTONIC, &end);

    mass_radius(s_gp, mu, &eosBM, &evolution_functions, r_ratio_BM, r_e_BM, Omega_BM, &Mass_BM, &Mass_0_BM, &J_BM, &R_e_BM, v_plus, v_minus, &Omega_K_BM, &star_props, &eosDM, r_ratio_DM, r_e_DM, Omega_DM, &Mass_DM, &Mass_0_DM, &J_DM, &R_e_DM, &DM_props, MDM, &R_p_BM, &R_p_DM, &Omega_K_DM, &T_BM, &T_DM, &W_BM, &W_DM, counter, &GRV2, &GRV3);

    fdm_old = fdm;
    fdm = Mass_DM / (Mass_BM + Mass_DM);

    if (verbose >= STANDARD) {
        double time_spent = (end.tv_sec - start.tv_sec) +
                            (end.tv_nsec - start.tv_nsec) / 1e9;
        printf("Time taken: %g seconds\n", time_spent);
        printf("New BM: A = %g   B = %g\n", DiffRotBM.A, DiffRotBM.B);
        printf("New DM: A = %g   B = %g\n", DiffRotDM.A, DiffRotDM.B);
        printf("\nGRV2 = %g  (Eq. 61 of Astron. Astrophys. Suppl. Ser. 132, 431-454) \n", GRV2);
        printf("GRV3 = %g  (Eq. 62 of Astron. Astrophys. Suppl. Ser. 132, 431-454) \n", GRV3);
        printf("Found Dark Matter fraction: %g\n\n", fdm);
    }

    if (verbose >= STANDARD) {
        printf("Baryonic Matter:\n");
        printf("  %-15s %-15s %-15s %-15s\n", "Mass (Msun)", "Mass_0 (Msun)", "R_e (Km)", "J / M^2");
        printf("  %-15.8g %-15.8g %-15.8f %-15.8f\n",
               Mass_BM / MSUN, Mass_0_BM / MSUN, R_e_BM / 1.0e5, (C * J_BM / (G * Mass_BM * Mass_BM)));
        if (!zero) {
            printf("\nDark Matter:\n");
            printf("  %-15s %-15s %-15s %-15s %-15s\n", "Mass (Msun)", "Mass_0 (Msun)", "R_e (Km)", "J / M^2", "fdm");
            printf("  %-15.8g %-15.8g %-15.8f %-15.8f %-15.8f\n",
                   Mass_DM / MSUN, Mass_0_DM / MSUN, R_e_DM / 1.0e5, (C * J_DM / (G * Mass_DM * Mass_DM + DBL_EPSILON)), fdm);
        }
        printf("\n");
    }

    FILE* file;
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

        // Extract EOS name
        char eos_name[128] = "unknownEOS";
        if (eosBM.type == TAB) {
            char* last_slash = strrchr(eosBM.filepath, '/');
            if (last_slash != NULL) {
                strncpy(eos_name, last_slash + 1, sizeof(eos_name));
                eos_name[sizeof(eos_name) - 1] = '\0';
                char* dot = strrchr(eos_name, '.');
                if (dot) *dot = '\0';  // Remove file extension if any
            }
        } else if (eosBM.type == POLY) {
            strcpy(eos_name, "POLY");
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