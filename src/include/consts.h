/***************************************************************************
 * consts.h
 *
 * This file contains definitions of constants used in the procedures.
 *
 **************************************************************************/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#define SLIM (2 * SDIV / 3)                  /* limit to speed up Omega computations */
#define DM (1.0 / (MDIV - 1.0))              /* spacing in mu direction */
#define RDIV 900                             /* grid point in RK integration */
#define SMAX 0.9999                          /* maximum value of s-coordinate */
#define DS (SMAX / (SDIV - 1.0))             /* spacing in s-direction */
#define LMAX 10                              /* max. term in Legendre poly. */
#define C 2.99792458e10                      /* speed of light in vacuum cm/s */
#define C_km (C / 100000)                    /* speed of light in vacuum km/s */
#define G 6.6730831e-8                       /* gravitational constant cm^3 / g / s^2 */
#define KAPPA (1.0e-15 * C * C / G)          /* scaling factor */
#define KSCALE (KAPPA * G / (C * C * C * C)) /* another scaling factor */
#define MSUN 1.988435e33                     /* Mass of Sun in g */
#define SQ(x) ((x) * (x))                    /* square macro */
#define MB 1.66e-24                          /* baryon mass */
#define RMIN 1.0e-15                         /* use approximate TOV equations when \
                                                computing spherical star and r<RMIN */
#define MeVc2_to_grams 1.78266269594644E-27  /* MeV/c^2 to grams */

#ifndef PI
#define PI 3.141592653589793 /* what else */
#endif

#define NTRIES 4 /* Maximum number of tries in root finding */
#define M_MAX 6   /* Maximum l=m mode considered */
#define BOUND 1   /* Excludes last grid point, since we set everything \
                     to zero there */
#define MAXIT 300 /* Maximum number of iterations in rtsec_G */
#include <float.h>
#include <limits.h>

#define GZIP_COMPRESSION_LEVEL 6

#define MAX_TAB_SIZE 16001
#define PR 0

#ifndef DBL_EPSILON
#define DBL_EPSILON 1e-15
#endif
#define EPS DBL_EPSILON

#ifdef DDEBUG
  #define DPRINT(...) fprintf(stderr, __VA_ARGS__)
#else
  #define DPRINT(...) ((void)0)
#endif

enum Verbosity {
    SILENT,
    STANDARD,
    VERBOSE,
    VERYVERBOSE,
    DEBUG
};

enum EOSType {
    TAB,
    POLY,
    UNKNOWN
};

enum DiffRotType {
    UNIFORM,
    DIFF,
    U8,
    U9EXT,
    UNKNOWN_ROT
};

struct DiffRotParams {
    double A;
    double B;
    double lambda1;
    double lambda2;
    double p;
    double csi;
};

struct EOS {
    char filepath[80];
    enum EOSType type;
    /* Tabulated EoSs */
    double log_e_tab[MAX_TAB_SIZE];
    double log_p_tab[MAX_TAB_SIZE];
    double log_h_tab[MAX_TAB_SIZE];
    double log_n0_tab[MAX_TAB_SIZE];
    int n_tab;
    /* Poly EoSs*/
    double n_P;
    double Gamma_P;
};

struct stellar_properties {
    /* All the values at the center and surface of the choosen Neutron Star */
    double e_center;
    double p_center;
    double h_center;
    double e_surface;
    double p_surface;
    double enthalpy_min;
    enum DiffRotType RotType;
};

struct evolution_variables {
    double **rho;
    double **gama;
    double **alpha;
    double **omega;
    double **energyBM;
    double **pressureBM;
    double **enthalpyBM;
    double **velocity_sqBM;
    double **Omega_diffBM;
    double **energyDM;
    double **pressureDM;
    double **enthalpyDM;
    double **velocity_sqDM;
    double **Omega_diffDM;
};

typedef void (*UpdateFunc)(const char *value, void *var);

// Structure to hold key, update function, and variable pointer
typedef struct
{
    const char *key;
    UpdateFunc update_func;
    void *variable;
} ConfigEntry;

struct rk_coeff {
    double a;
    double b;
    double b_DM;
    double c;
    double c_DM;
};

#endif