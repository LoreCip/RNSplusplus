#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "consts.h"
#include "equil.h"
#include "equil_util.h"
#include "nrutil.h"
#include "output.h"
#include "rnsregrid_util.h"

#define H5_USE_18_API
#define H5Acreate_vers 2
#define H5Dcreate_vers 2
#define H5Dopen_vers 2
#include "hdf5.h"

int main(int argc, char **argv) {
    struct EOS eosBM;
    struct EOS eosDM;
    struct stellar_properties star_props;
    struct stellar_properties DM_props;

    struct evolution_variables ev;
    allocate_evolution_variables(&ev);

    double r_ratio_BM, /* axis ratio */
        r_ratio_DM,
        Omega_BM,
        Omega_DM,
        r_e_BM,
        r_e_DM,
        s_gp_tmp[SDIV + 1], s_gp[SDIV + 1], /* s grid points */
        mu_tmp[SDIV + 1], mu[MDIV + 1];     /* \mu grid points */

    double Mass_BM, Mass_DM, Mass_0_BM, Mass_0_DM, J_BM, J_DM, R_e_BM, R_e_DM, v_plus, v_minus, Omega_K;

    /* ******************************************
       Default Values
       ********************************************* */

    double Omega_pt[1], R_e_pt[1], r_e_pt[1], mass0_pt[1];
    double rho0_center = 1.28e-3; /* Central density */
    char eos_type[80] = "poly";   /* {"poly" "tab"} */
    char eos_file[80] = "no EOS file specified";
    double eos_k = 100.0;               /* poly EOS K */
    double Gamma_P = 2.0;               /* poly EOS Gamma */
    char rotation_type[20] = "uniform"; /* {"uniform"} */
    double r_ratio = 1.00;              /* axes ratio r_p/r_e */
    double accuracy = 1.e-7;            /* accuracy goal */
    double atm_factor = 1.e-10;

    char zero_shift[20] = "no";
    char save_2Dmodel[20] = "no";
    char model2D_file_IN[256] = "input_model.h5";
    char model_file_OUT[256] = "output.h5";
    char recover_2Dmodel[20] = "no";
    double cf = 1.0;
    int RNS_lmax = 10;

    // Delfault diff. rot params
    struct DiffRotParams DiffRotBM, DiffRotDM;
    DiffRotBM.A = 0.8;
    DiffRotBM.B = 0.6;
    DiffRotBM.lambda1 = 2.0;
    DiffRotBM.lambda2 = 0.5;
    DiffRotBM.p = 1;
    DiffRotBM.csi = 3;

    DiffRotDM.A = 0.8;
    DiffRotDM.B = 0.6;
    DiffRotDM.lambda1 = 2.0;
    DiffRotDM.lambda2 = 0.5;
    DiffRotDM.p = 1;
    DiffRotDM.csi = 3;


    // -----------------------------
    //  Carteisan GRID size  (default values)
    // -----------------------------
    int np;  // Total number of grid points
    int Px, Py, Pz, Nx, Ny, Nz;
    double Lx, Ly, Lz, dx;
    Px = Py = Pz = 25;
    Nx = Ny = Nz = 49;
    Lx = Ly = Lz = 12;
    dx = 0.5;
    np = Nx * Ny * Nz;

    /****************************************************************** */

    int opterr, c;
    opterr = 0;
    for (c = getopt(argc, argv, "+ht:f:k:g:i:o:x:y:z:d:");
         c != EOF;
         c = getopt(argc, argv, "+ht:f:k:g:i:o:x:y:z:d:")) {
        switch (c) {
            case 'x': /* Number of grid point x-direction from 0 (included) to BORDER */
                sscanf(optarg, "%d", &Px);
                break;
            case 'y': /* Number of grid point x-direction from 0 (included) to BORDER */
                sscanf(optarg, "%d", &Py);
                break;
            case 'z': /* Number of grid point x-direction from 0 (included) to BORDER */
                sscanf(optarg, "%d", &Pz);
                break;
            case 'd': /* grid spacing in units G=c=M_sun=1 */
                sscanf(optarg, "%lf", &dx);
                break;
            case 'i': /* initial_model_data.h5  */
                sscanf(optarg, "%s", model2D_file_IN);
                break;
            case 'o': /* ginal_model_data.h5 */
                sscanf(optarg, "%s", model_file_OUT);
                break;
            case 't': /* eos type */
                sscanf(optarg, "%s", eos_type);
                break;
            case 'f': /* eos TAB file */
                sscanf(optarg, "%s", eos_file);
                break;
            case 'k': /* poly EOS K */
                sscanf(optarg, "%lf", &eos_k);
                break;
            case 'g': /* poly EOS Gamma */
                sscanf(optarg, "%lf", &Gamma_P);
                break;
            case 'h':
                printf("\n");
                printf("Help: ./RNS_readID -options<meaning> {default value} \n");
                printf("\n");
                printf("      -i <initial_model_data.hdf> \n");
                printf("      -o <initial_model_grid_data.hdf> \n");
                printf("\n");
                printf("      -d <step-in-coordinates> \n");
                printf("      -x <number of points in x-dir>  {12} \n");
                printf("      -y <number of points in x-dir>  {12} \n");
                printf("      -z <number of points in x-dir>  {12} \n");
                printf("      [GRID size, #of points = (2*Px-1) * (2*Py-1) * (2*Pz-1)]\n");
                printf("      [coordinates = d [-Pz,Pz]) x (d [-Py,Py])x (d [-Px,Px])]\n");
                printf("\n");
                return 0;
            default:
                break;
        }
    }

    /*********************************************************** */

    int sdiv;
    int mdiv;
    printf("Reading data from file %s\n", model2D_file_IN);

    hdf5_read_var(s_gp_tmp, &sdiv, mu_tmp, &mdiv, model2D_file_IN, &ev, &star_props, &DM_props, &eosBM, &eosDM, &r_ratio_BM, &r_ratio_DM, &Omega_BM, &Omega_DM, &r_e_BM, &r_e_DM, &DiffRotBM, &DiffRotDM);
    for (int s = 1; s <= SDIV; s++){
        s_gp[s] = s_gp_tmp[s - 1];
    }
    for (int m = 1; m <= MDIV; m++){
        mu[m] = mu_tmp[m - 1];
    }

    /*
        Define Cartesian Grid
    */
    Nx = 2 * Px - 1;
    Ny = 2 * Py - 1;
    Nz = 2 * Pz - 1;
    Lx = (Px - 1) * dx;
    Ly = (Py - 1) * dx;
    Lz = (Pz - 1) * dx;
    np = Nx * Ny * Nz;
    printf("\n====================================================\n");
    printf(" Grid setting x range is: [ %e , %e ]\n", -Lx, Lx);
    printf("              y range is: [ %e , %e ]\n", -Ly, Ly);
    printf("              z range is: [ %e , %e ]\n", -Lz, Lz);
    printf(" Number of points is (%d) %d,%d,%d \n", np, Nz, Ny, Nx);
    printf("====================================================\n");

    //**************************************************************************************************
    /* CONSTRUCT CARTESIAN ARRAYS WITH NEEDED POLAR QUANTITIES */
    //**************************************************************************************************

    /// 1-D array storing the values of coordinate x of the {\tt np} grid points [unit: km]
    double *x_grid;
    /// 1-D array storing the values of coordinate y of the {\tt np} grid points [unit: km]
    double *y_grid;
    /// 1-D array storing the values of coordinate z of the {\tt np} grid points [unit: km]
    double *z_grid;
    /// Lapse function $N$ at the {\tt np} grid points (1-D array)
    double *alp;    /// Component $\beta^x$ of the shift vector of non rotating coordinates [unit: $c$]
    double *betax;  /// Component $\beta^y$ of the shift vector of non rotating coordinates [unit: $c$]
    double *betay;  /// Component $\beta^z$ of the shift vector of non rotating coordinates [unit: $c$]
    double *betaz;  /// Metric coefficient $\gamma_{xx}$ at the grid points (1-D array)
    double *gxx;
    /// Metric coefficient $\gamma_{xy}$ at the grid points (1-D array)
    double *gxy;
    /// Metric coefficient $\gamma_{xz}$ at the grid points (1-D array)
    double *gxz;
    /// Metric coefficient $\gamma_{yy}$ at the grid points (1-D array)
    double *gyy;
    /// Metric coefficient $\gamma_{yz}$ at the grid points (1-D array)
    double *gyz;
    /// Metric coefficient $\gamma_{zz}$ at the grid points (1-D array)
    double *gzz;
    /// Component $K_{xx}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
    double *kxx;
    /// Component $K_{xy}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
    double *kxy;
    /// Component $K_{xz}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
    double *kxz;
    /// Component $K_{yy}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
    double *kyy;
    /// Component $K_{yz}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
    double *kyz;
    /// Component $K_{zz}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
    double *kzz;
    // Hydro components
    //------------------
    /** Baryon density in the fluid frame at the {\tt np} grid points (1-D array)
     * [unit: ${\rm kg \, m}^{-3}$]
     */
    double *rho_BM, *rho_DM;
    double *press_BM, *press_DM;
    /// Specific internal energy at the  {\tt np} grid points (1-D array) [unit: $c^2$]
    double *eps_BM, *eps_DM;
    /** Component $U^x$ of the fluid 3-velocity with respect to the Eulerian
     * observer, at the {\tt np} grid points (1-D array) [unit: $c$]
     */
    double *velx_BM, *velx_DM;
    /** Component $U^y$ of the fluid 3-velocity with respect to the Eulerian
     * observer, at the {\tt np} grid points (1-D array) [unit: $c$]
     */
    double *vely_BM, *vely_DM;
    /** Component $U^z$ of the fluid 3-velocity with respect to the Eulerian
     * observer, at the {\tt np} grid points (1-D array) [unit: $c$]
     */
    double *velz_BM, *velz_DM;
    double *Odif_BM, *Odif_DM;

    int ii, jj, kk, idx, s, m;
    double
        n_P,              /* polytropic index */
        **nu,             /* potential nu */
        **B,              /* potential B */
        **rho_0_BM,       /* rest mass density */
        **rho_0_DM,       /* rest mass density */
        **nu_dr,          /* r-der. in s-coord. of nu */
        **B_dr,           /* r-der. in s-coord. of B */
        **alpha_dr,       /* r-der. in s-coord. of alpha */
        **omega_dr,       /* r-der. in s-coord. of omega */
        **nu_dth,         /* theta-der. in mu-coord. of nu */
        **B_dth,          /* theta-der. in mu-coord. of B */
        **alpha_dth,      /* theta-der. in mu-coord. of alpha */
        **omega_dth,      /* theta-der. in mu-coord. of omega */
        x_i,              /* x at i */
        y_j,              /* y at j */
        z_k,              /* z at k */
        nu_ijk,           /* nu at ijk point */
        exp_nu_ijk,       /* exp(nu) at ijk point */
        B_ijk,            /* B at ijk point */
        omega_ijk,        /* omega at ijk point */
        alpha_ijk,        /* alpha at ijk point */
        exp_alpha_ijk,    /* exp(alpha) at ijk point */
        rho_0_BM_ijk,     /* rho_0 at ijk point */
        energy_BM_ijk,    /* energy at ijk point */
        pressure_BM_ijk,  /* pressure at ijk point */
        rho_0_DM_ijk,     /* rho_0 at ijk point */
        energy_DM_ijk,    /* energy at ijk point */
        pressure_DM_ijk,  /* pressure at ijk point */
        nu_dx,            /* derivative of nu w.r.t. x */
        nu_dy,            /* derivative of nu w.r.t. y */
        B_dx,             /* derivative of B w.r.t. x */
        B_dy,             /* derivative of B w.r.t. y */
        omega_dx,         /* derivative of omega w.r.t. x */
        omega_dy,         /* derivative of omega w.r.t. y */
        omega_dz,         /* derivative of omega w.r.t. z */
        alpha_dx,         /* derivative of alpha w.r.t. x */
        alpha_dy,         /* derivative of alpha w.r.t. y */
        r_ijk,            /* r at ijk point */
        r_bar_ijk,        /* sqrt(x^2+y^2) at ijk point */
        dr_dx,            /* dr/dx */
        dr_dy,            /* dr/dy */
        dr_dz,            /* dr/dz */
        dtheta_dx,        /* dtheta/dx */
        dtheta_dy,        /* dtheta/dy */
        dtheta_dz,        /* dtheta/dz */
        nu_dr_ijk,        /* dnu/dr at ijk */
        B_dr_ijk,         /* dB/dr at ijk */
        alpha_dr_ijk,     /* dalpha/dr at ijk */
        omega_dr_ijk,     /* domega/dr at ijk */
        nu_dtheta_ijk,    /* dnu/dtheta at ijk */
        B_dtheta_ijk,     /* dB/dtheta at ijk */
        alpha_dtheta_ijk, /* dalpha/dtheta at ijk */
        omega_dtheta_ijk, /* domega/dtheta at ijk */
        gamma_ijk,        /* gamma = det(3g) */
        h_BM_ijk,         /* h = 1 + eps + P/rho */
        h_DM_ijk,         /* h = 1 + eps + P/rho */
        W_BM_ijk,         /* Lorentz factor */
        W_DM_ijk,         /* Lorentz factor */
        distance_ijk,     /* Signed distance to surface */
        e_atm,            /* energy density of atmosphere */
        p_atm,            /* pressure of atmosphere */
        rho_0_atm,        /* rest mass density of atmosphere */
        dens_atm,         /* D of atmosphere */
        tau_atm,          /* tau of atmosphere */
        temp_a,           /* temporary variables */
        temp_o,
        temp_g,
        temp_r,
        temp_e,
        temp_p,
        temp_h,
        temp_v,
        OmegaBM_ijk,
        OmegaDM_ijk;

    x_grid = malloc(np * sizeof(double));  // new double[np] ;
    y_grid = malloc(np * sizeof(double));  // new double[np] ;
    z_grid = malloc(np * sizeof(double));  // new double[np] ;

    alp = malloc(np * sizeof(double));  // new double[np] ;

    betax = malloc(np * sizeof(double));  // new double[np] ;
    betay = malloc(np * sizeof(double));  // new double[np] ;
    betaz = malloc(np * sizeof(double));  // new double[np] ;

    gxx = malloc(np * sizeof(double));  // new double[np] ;
    gxy = malloc(np * sizeof(double));  // new double[np] ;
    gxz = malloc(np * sizeof(double));  // new double[np] ;
    gyy = malloc(np * sizeof(double));  // new double[np] ;
    gyz = malloc(np * sizeof(double));  // new double[np] ;
    gzz = malloc(np * sizeof(double));  // new double[np] ;

    kxx = malloc(np * sizeof(double));  // new double[np] ;
    kxy = malloc(np * sizeof(double));  // new double[np] ;
    kxz = malloc(np * sizeof(double));  // new double[np] ;
    kyy = malloc(np * sizeof(double));  // new double[np] ;
    kyz = malloc(np * sizeof(double));  // new double[np] ;
    kzz = malloc(np * sizeof(double));  // new double[np] ;

    rho_BM = malloc(np * sizeof(double));    // new double[np] ;
    press_BM = malloc(np * sizeof(double));  // new double[np] ;
    eps_BM = malloc(np * sizeof(double));    // new double[np] ;

    velx_BM = malloc(np * sizeof(double));  // new double[np] ;
    vely_BM = malloc(np * sizeof(double));  // new double[np] ;
    velz_BM = malloc(np * sizeof(double));  // new double[np] ;

    Odif_BM = malloc(np * sizeof(double));  // new double[np] ;

    rho_DM = malloc(np * sizeof(double));    // new double[np] ;
    press_DM = malloc(np * sizeof(double));  // new double[np] ;
    eps_DM = malloc(np * sizeof(double));    // new double[np] ;

    velx_DM = malloc(np * sizeof(double));  // new double[np] ;
    vely_DM = malloc(np * sizeof(double));  // new double[np] ;
    velz_DM = malloc(np * sizeof(double));  // new double[np] ;

    Odif_DM = malloc(np * sizeof(double));  // new double[np] ;

    double r_e = FMAX(r_e_BM, r_e_DM) * sqrt(KAPPA) / 1e5;

    for (ii = 0; ii < Nx; ii++)
        for (jj = 0; jj < Ny; jj++)
            for (kk = 0; kk < Nz; kk++) {
                idx = ii + jj * Nx + kk * Nx * Ny;
                x_grid[idx] = -Lx + dx * ii;
                y_grid[idx] = -Ly + dx * jj;
                z_grid[idx] = -Lz + dx * kk;
            }

    //**************************************************************************************************
    //**************************************************************************************************

    nu = array_allocate(1, SDIV, 1, MDIV);
    B = array_allocate(1, SDIV, 1, MDIV);
    rho_0_BM = array_allocate(1, SDIV, 1, MDIV);
    rho_0_DM = array_allocate(1, SDIV, 1, MDIV);

    for (m = 1; m <= MDIV; m++) {
        for (s = 1; s <= SDIV; s++) {
            nu[s][m] = (ev.gama[s][m] + ev.rho[s][m]) / 2.0;
            B[s][m] = exp(ev.gama[s][m]);
            rho_0_BM[s][m] = (ev.energyBM[s][m] + ev.pressureBM[s][m]) * exp(-ev.enthalpyBM[s][m]);
            rho_0_DM[s][m] = (ev.energyDM[s][m] + ev.pressureDM[s][m]) * exp(-ev.enthalpyDM[s][m]);
        }
    }

    array_free(ev.rho, 1, SDIV, 1, MDIV);
    array_free(ev.gama, 1, SDIV, 1, MDIV);
    array_free(ev.enthalpyBM, 1, SDIV, 1, MDIV);
    array_free(ev.velocity_sqBM, 1, SDIV, 1, MDIV);
    array_free(ev.enthalpyDM, 1, SDIV, 1, MDIV);
    array_free(ev.velocity_sqDM, 1, SDIV, 1, MDIV);

    nu_dr = array_allocate(1, SDIV, 1, MDIV);
    B_dr = array_allocate(1, SDIV, 1, MDIV);
    alpha_dr = array_allocate(1, SDIV, 1, MDIV);
    omega_dr = array_allocate(1, SDIV, 1, MDIV);
    nu_dth = array_allocate(1, SDIV, 1, MDIV);
    B_dth = array_allocate(1, SDIV, 1, MDIV);
    alpha_dth = array_allocate(1, SDIV, 1, MDIV);
    omega_dth = array_allocate(1, SDIV, 1, MDIV);

    for (s = 1; s <= SDIV; s++) {
        for (m = 1; m <= MDIV; m++) {
            nu_dr[s][m] = deriv_s(nu, s, m) * SQ(1.0 - s_gp[s]) / r_e;
            B_dr[s][m] = deriv_s(B, s, m) * SQ(1.0 - s_gp[s]) / r_e;
            alpha_dr[s][m] = deriv_s(ev.alpha, s, m) * SQ(1.0 - s_gp[s]) / r_e;
            omega_dr[s][m] = deriv_s(ev.omega, s, m) * SQ(1.0 - s_gp[s]) / r_e;
            nu_dth[s][m] = deriv_m(nu, s, m) * (-sqrt(1.0 - SQ(mu[m])));
            B_dth[s][m] = deriv_m(B, s, m) * (-sqrt(1.0 - SQ(mu[m])));
            alpha_dth[s][m] = deriv_m(ev.alpha, s, m) * (-sqrt(1.0 - SQ(mu[m])));
            omega_dth[s][m] = deriv_m(ev.omega, s, m) * (-sqrt(1.0 - SQ(mu[m])));
        }
    }

    /* COMPUTE INITIAL DATA */

    rho_0_atm = 1e-10; /* rename the constant for historical reasons */
    e_atm = rho_0_atm;
    p_atm = eos_k * pow(rho_0_atm, Gamma_P);

    printf("Computing cartesian grids...\n");
    fflush(stdout);

    int i, j, k;
// #pragma omp parallel for collapse(3) private(i, j, k, idx, x_i, y_j, z_k, nu_ijk, exp_nu_ijk, B_ijk, omega_ijk, alpha_ijk, exp_alpha_ijk, rho_0_BM_ijk, energy_BM_ijk, pressure_BM_ijk, rho_0_DM_ijk, energy_DM_ijk, pressure_DM_ijk, OmegaBM_ijk, OmegaDM_ijk, r_ijk, r_bar_ijk, dr_dx, dr_dy, dr_dz, dtheta_dx, dtheta_dy, dtheta_dz, nu_dx, nu_dy, B_dx, B_dy, alpha_dx, alpha_dy, omega_dx, omega_dy, omega_dz, gamma_ijk, h_BM_ijk, h_DM_ijk, W_BM_ijk, W_DM_ijk, distance_ijk, dens_atm, tau_atm)
    for (i = 1; i <= Nx; i++) {
        for (j = 1; j <= Ny; j++) {
            for (k = 1; k <= Nz; k++) {
                int idx = i - 1 + Nx * (j - 1 + Ny * (k - 1));
                x_i = x_grid[idx];
                y_j = y_grid[idx];
                z_k = z_grid[idx];

                grid_interp_all(s_gp, mu, r_e, Nx, Ny, Nz,
                                x_grid, y_grid, z_grid,
                                i, j, k,
                                nu, B, ev.alpha, ev.omega,
                                nu_dr, B_dr, alpha_dr, omega_dr,
                                nu_dth, B_dth, alpha_dth, omega_dth,
                                rho_0_BM, ev.energyBM, ev.pressureBM,
                                rho_0_DM, ev.energyDM, ev.pressureDM,
                                ev.Omega_diffBM, ev.Omega_diffDM,
                                &nu_ijk, &B_ijk, &alpha_ijk, &omega_ijk,
                                &nu_dr_ijk, &B_dr_ijk, &alpha_dr_ijk,
                                &omega_dr_ijk,
                                &nu_dtheta_ijk, &B_dtheta_ijk,
                                &alpha_dtheta_ijk, &omega_dtheta_ijk,
                                &rho_0_BM_ijk, &energy_BM_ijk, &pressure_BM_ijk,
                                &rho_0_DM_ijk, &energy_DM_ijk, &pressure_DM_ijk,
                                &distance_ijk, &OmegaBM_ijk, &OmegaDM_ijk
                );

                /* *************************************** */
                /* DETECT if it is in the ATMOSPHERE       */
                /* *************************************** */
                if ((rho_0_BM_ijk <= 0.0) || (energy_BM_ijk <= 0.0) ||
                    (pressure_BM_ijk <= 0.0)) {
                    rho_0_BM_ijk = rho_0_atm;
                    energy_BM_ijk = e_atm + 1.e-20;
                    pressure_BM_ijk = p_atm;
                }
                if ((rho_0_DM_ijk <= 0.0) || (energy_DM_ijk <= 0.0) ||
                    (pressure_DM_ijk <= 0.0)) {
                    rho_0_DM_ijk = rho_0_atm;
                    energy_DM_ijk = e_atm + 1.e-20;
                    pressure_DM_ijk = p_atm;
                }
                /* *************************************** */
                /* END ATMOSPERE SETTINGs                      */
                /* *************************************** */

                exp_nu_ijk = exp(nu_ijk);
                exp_alpha_ijk = exp(alpha_ijk);

                r_ijk = sqrt(SQ(x_i) + SQ(y_j) + SQ(z_k));
                r_bar_ijk = sqrt(SQ(x_i) + SQ(y_j));

                alp[idx] = exp_nu_ijk;

                if (x_i == 0.0 && y_j == 0.0) {
                    gxx[idx] = SQ(exp_alpha_ijk);
                    gyy[idx] = SQ(exp_alpha_ijk);
                    gzz[idx] = SQ(exp_alpha_ijk);

                    gxy[idx] = 0.0;
                    gxz[idx] = 0.0;
                    gyz[idx] = 0.0;

                    kxx[idx] = 0.0;
                    kyy[idx] = 0.0;
                    kzz[idx] = 0.0;
                    kxy[idx] = 0.0;
                    kxz[idx] = 0.0;
                    kyz[idx] = 0.0;
                } else {
                    dr_dx = x_i / r_ijk;
                    dr_dy = y_j / r_ijk;
                    dr_dz = z_k / r_ijk;

                    dtheta_dx = x_i * z_k / (SQ(r_ijk) * r_bar_ijk);
                    dtheta_dy = y_j * z_k / (SQ(r_ijk) * r_bar_ijk);
                    dtheta_dz = -r_bar_ijk / SQ(r_ijk);

                    nu_dx = dr_dx * nu_dr_ijk + dtheta_dx * nu_dtheta_ijk;
                    nu_dy = dr_dy * nu_dr_ijk + dtheta_dy * nu_dtheta_ijk;

                    B_dx = dr_dx * B_dr_ijk + dtheta_dx * B_dtheta_ijk;
                    B_dy = dr_dy * B_dr_ijk + dtheta_dy * B_dtheta_ijk;

                    alpha_dx = dr_dx * alpha_dr_ijk + dtheta_dx * alpha_dtheta_ijk;
                    alpha_dy = dr_dy * alpha_dr_ijk + dtheta_dy * alpha_dtheta_ijk;

                    omega_dx = dr_dx * omega_dr_ijk + dtheta_dx * omega_dtheta_ijk;
                    omega_dy = dr_dy * omega_dr_ijk + dtheta_dy * omega_dtheta_ijk;

                    /* enforce omega_dz=0 at z=0 (it is slightly nonzero due
                       to O(h) forwards formula in computing derivative) */

                    if (z_k == 0.0)
                        omega_dz = 0.0;
                    else
                        omega_dz = dr_dz * omega_dr_ijk + dtheta_dz * omega_dtheta_ijk;

                    gxx[idx] = (SQ(B_ijk * y_j / exp_nu_ijk) + SQ(exp_alpha_ijk * x_i)) /
                               (SQ(x_i) + SQ(y_j));

                    gxy[idx] = (SQ(exp_alpha_ijk) - SQ(B_ijk / exp_nu_ijk)) *
                               x_i * y_j / (SQ(x_i) + SQ(y_j));

                    gxz[idx] = 0.0;

                    gyy[idx] = (SQ(B_ijk * x_i / exp_nu_ijk) + SQ(exp_alpha_ijk * y_j)) /
                               (SQ(x_i) + SQ(y_j));

                    gyz[idx] = 0.0;

                    gzz[idx] = SQ(exp_alpha_ijk);

                    kxx[idx] = ((SQ(r_bar_ijk) * y_j * omega_dx +
                                 (x_i * nu_dy - y_j * nu_dx) * SQ(y_j) * omega_ijk) *
                                    SQ(B_ijk) +
                                (y_j * B_dx - x_i * B_dy) * omega_ijk * SQ(y_j) * B_ijk + (y_j * alpha_dx - x_i * alpha_dy) * omega_ijk * SQ(x_i * exp_alpha_ijk * exp_nu_ijk)) /
                               (SQ(r_bar_ijk * exp_nu_ijk) * exp_nu_ijk);

                    kxy[idx] = ((0.5 * SQ(r_bar_ijk) *
                                     (y_j * omega_dy - x_i * omega_dx) +
                                 (y_j * nu_dx - x_i * nu_dy) * x_i * y_j *
                                     omega_ijk) *
                                    SQ(B_ijk) +
                                (-y_j * B_dx + x_i * B_dy) * omega_ijk * x_i * y_j * B_ijk + (y_j * alpha_dx - x_i * alpha_dy) * omega_ijk * x_i * y_j * SQ(exp_alpha_ijk * exp_nu_ijk)) /
                               (SQ(r_bar_ijk * exp_nu_ijk) * exp_nu_ijk);

                    kxz[idx] = 0.5 * SQ(B_ijk) * y_j * omega_dz /
                               (SQ(exp_nu_ijk) * exp_nu_ijk);

                    kyy[idx] = ((-SQ(r_bar_ijk) * x_i * omega_dy +
                                 (x_i * nu_dy - y_j * nu_dx) * SQ(x_i) *
                                     omega_ijk) *
                                    SQ(B_ijk) +
                                (y_j * B_dx - x_i * B_dy) * omega_ijk * SQ(x_i) * B_ijk + (y_j * alpha_dx - x_i * alpha_dy) * omega_ijk * SQ(y_j * exp_alpha_ijk * exp_nu_ijk)) /
                               (SQ(r_bar_ijk * exp_nu_ijk) * exp_nu_ijk);

                    kyz[idx] = -0.5 * SQ(B_ijk) * x_i * omega_dz /
                               (SQ(exp_nu_ijk) * exp_nu_ijk);

                    kzz[idx] = (y_j * alpha_dx - x_i * alpha_dy) *
                               omega_ijk * SQ(exp_alpha_ijk) /
                               exp_nu_ijk;
                }

                betax[idx] = omega_ijk * y_j;
                betay[idx] = -omega_ijk * x_i;

                betaz[idx] = 0.0;

                rho_BM[idx] = rho_0_BM_ijk;
                rho_DM[idx] = rho_0_DM_ijk;

                eps_BM[idx] = energy_BM_ijk / rho_0_BM_ijk - 1.0;
                eps_DM[idx] = energy_DM_ijk / rho_0_DM_ijk - 1.0;

                h_BM_ijk = (energy_BM_ijk + pressure_BM_ijk) / rho_0_BM_ijk;
                h_DM_ijk = (energy_DM_ijk + pressure_DM_ijk) / rho_0_DM_ijk;

                gamma_ijk = SQ(exp_alpha_ijk) * SQ(B_ijk * exp_alpha_ijk /
                                                   exp_nu_ijk);

                W_BM_ijk = 1.0 / sqrt(1.0 - SQ((omega_ijk - OmegaBM_ijk) * B_ijk *
                                               r_bar_ijk / SQ(exp_nu_ijk)));
                W_DM_ijk = 1.0 / sqrt(1.0 - SQ((omega_ijk - OmegaDM_ijk) * B_ijk *
                                               r_bar_ijk / SQ(exp_nu_ijk)));

                velx_BM[idx] = (omega_ijk - OmegaBM_ijk) * y_j / exp_nu_ijk / C_km;
                velx_DM[idx] = (omega_ijk - OmegaDM_ijk) * y_j / exp_nu_ijk / C_km;

                vely_BM[idx] = -(omega_ijk - OmegaBM_ijk) * x_i / exp_nu_ijk / C_km;
                vely_DM[idx] = -(omega_ijk - OmegaDM_ijk) * x_i / exp_nu_ijk / C_km;

                velz_BM[idx] = 0.0;
                velz_DM[idx] = 0.0;

                press_BM[idx] = pressure_BM_ijk;
                press_DM[idx] = pressure_DM_ijk;

                Odif_BM[idx] = OmegaBM_ijk;
                Odif_DM[idx] = OmegaDM_ijk;

                dens_atm = sqrt(gamma_ijk) * rho_0_atm;
                tau_atm = sqrt(gamma_ijk) * eos_k * pow(rho_0_atm, Gamma_P) /
                          (Gamma_P - 1.0);

                /* *************************************** */
                /* ATMOSPERE SETTINGs                      */
                /* *************************************** */
                if ((rho_BM[idx] < (1.0 + 1e-3) * rho_0_atm)) {
                    velx_BM[idx] = 0.0;
                    vely_BM[idx] = 0.0;
                    velz_BM[idx] = 0.0;
                }
                if ((rho_DM[idx] < (1.0 + 1e-3) * rho_0_atm)) {
                    velx_DM[idx] = 0.0;
                    vely_DM[idx] = 0.0;
                    velz_DM[idx] = 0.0;
                }
                /* *************************************** */
                /* ATMOSPERE SETTINGs                      */
                /* *************************************** */
            }
        }
    }

    // ------------------------------------------
    // HDF5 save
    // ------------------------------------------
    printf("Creating file. This can be very slow!\n");
    fflush(stdout);

    hid_t file_id; /* file identifier */
    herr_t status;

    hid_t dataset_id, dataspace_id, plist_id; /* identifiers for dsets*/
    hsize_t dims[3], chunk_dims[3];

    double *dset_data;
    double **var;
    int varIndex;

    /* =========================================== */
    /* Create a new file using default properties. */
    /* =========================================== */
    file_id = H5Fcreate(model_file_OUT, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dims[2] = Nx;
    dims[1] = Ny;
    dims[0] = Nz;
    calculate_chunk_dims(chunk_dims, dims);

    double *varlist[] = {x_grid, y_grid, z_grid,
                         alp, betax, betay, betaz,
                         gxx, gxy, gxz, gyy, gyz, gzz,
                         kxx, kxy, kxz, kyy, kyz, kzz,
                         rho_BM, press_BM, eps_BM, velx_BM, vely_BM, velz_BM, Odif_BM,
                         rho_DM, press_DM, eps_DM, velx_DM, vely_DM, velz_DM, Odif_DM};
    char *varnames[] = {"/X[0]", "/X[1]", "/X[2]",
                        "/alp", "/betax", "/betay", "/betaz",
                        "/gxx", "/gxy", "/gxz", "/gyy", "/gyz", "/gzz",
                        "/Kxx", "/Kxy", "/Kxz", "/Kyy", "/Kyz", "/Kzz",
                        "/rhoBM", "/pressBM", "/epsBM",
                        "/velBM[0]", "/velBM[1]", "/velBM[2]", "/OmegaBM",
                        "/rhoDM", "/pressDM", "/epsDM",
                        "/velDM[0]", "/velDM[1]", "/velDM[2]", "/OmegaDM"};

    for (idx = 0; idx < 33; idx++) {
        printf("\rCurrent variable: %s (%d/33)%-5s", varnames[idx], idx + 1, "");
        fflush(stdout);

        dataspace_id = H5Screate_simple(3, dims, NULL);
        plist_id = H5Pcreate(H5P_DATASET_CREATE);

        H5Pset_chunk(plist_id, 3, chunk_dims);  // Set chunking
        H5Pset_deflate(plist_id, 6);            // Set gzip compression
        H5Pset_shuffle(plist_id);

        dataset_id = H5Dcreate(file_id, varnames[idx],
                               H5T_NATIVE_DOUBLE, dataspace_id,
                               H5P_DEFAULT, plist_id, H5P_DEFAULT);

        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          varlist[idx]);

        H5Pclose(plist_id);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
    }

    status = H5Fclose(file_id);

    // ------------------------------------------
    // Deallocate variables
    // ------------------------------------------

    array_free(ev.alpha, 1, SDIV, 1, MDIV);
    array_free(ev.omega, 1, SDIV, 1, MDIV);

    array_free(rho_0_BM, 1, SDIV, 1, MDIV);
    array_free(rho_0_DM, 1, SDIV, 1, MDIV);
    array_free(ev.energyBM, 1, SDIV, 1, MDIV);
    array_free(ev.pressureBM, 1, SDIV, 1, MDIV);
    array_free(ev.energyDM, 1, SDIV, 1, MDIV);
    array_free(ev.pressureDM, 1, SDIV, 1, MDIV);

    array_free(nu, 1, SDIV, 1, MDIV);
    array_free(B, 1, SDIV, 1, MDIV);

    array_free(nu_dr, 1, SDIV, 1, MDIV);
    array_free(B_dr, 1, SDIV, 1, MDIV);
    array_free(alpha_dr, 1, SDIV, 1, MDIV);
    array_free(omega_dr, 1, SDIV, 1, MDIV);
    array_free(nu_dth, 1, SDIV, 1, MDIV);
    array_free(B_dth, 1, SDIV, 1, MDIV);
    array_free(alpha_dth, 1, SDIV, 1, MDIV);
    array_free(omega_dth, 1, SDIV, 1, MDIV);

    printf("\nDone.\n");
}