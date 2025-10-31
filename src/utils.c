/**
 * utils.c
 *
 * Utility functions for galaxy/cluster simulations: Cosmology constants,
 * safe math operations, array allocation, unit conversions, and error handling.
 * Supports modules like density/profiles/sampling with astrophysical units.
 *
 * Cosmology: Planck18 approx (h=0.674, Om=0.315); extend for varying z.
 * Safe math: Wrappers for pow/log to avoid nan/inf/overflow.
 * Arrays: Dynamic allocation with error checks.
 *
 * Units: Primary: kpc, Msun, km/s; conversions to SI.
 *
 * References:
 * - Cosmology: Planck Collaboration (2020) [web:368][web:369]
 * - Math utils: Numerical Recipes (Press et al. 2007) [web:374]
 * - Constants: Colossus/Astropy [web:387][web:384]
 *
 * Usage: #include "utils.h"; link with -lm.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include "utils.h"  // Assumed header with prototypes

/* Cosmology constants (Planck18 defaults; editable) */
static double H0 = 67.4;          /* km/s/Mpc */
static double OMEGA_M = 0.315;    /* Matter density */
static double OMEGA_L = 0.685;    /* Lambda */
static double OMEGA_B = 0.049;    /* Baryon */
static double H = 0.674;          /* h = H0/100 */

/**
 * Compute critical density rho_crit(z=0) = 3 H0^2 / (8 pi G) [Msun / kpc^3]
 * @param h Optional h (default 0.674)
 * @return rho_crit [Msun kpc^{-3}]
 */
double compute_rho_crit(double h) {
    if (h <= 0.0) h = H;
    double H0_kmsMpc = H0 / h * h * 100.0;  /* Full H0 */
    double G_SI = 6.67430e-11;  /* m^3 kg^{-1} s^{-2} */
    double Msun_SI = 1.9885e30; /* kg */
    double kpc_m = 3.08568e19;  /* m */
    double rho_crit_SI = 3.0 * H0_kmsMpc * H0_kmsMpc * 1e5 * 1e5 / (8.0 * PI * G_SI);  /* kg m^{-3} */
    double rho_crit_ast = rho_crit_SI / (Msun_SI / pow(kpc_m, 3.0));  /* Msun kpc^{-3} */
    return rho_crit_ast * 1e-9;  /* Adjust scaling for h^2 */
}

/**
 * Safe power: pow(base, exp) with overflow/nan checks
 * @param base Base value
 * @param exp Exponent
 * @return result or 0.0 on error
 */
double safe_pow(double base, double exp) {
    if (isnan(base) || isnan(exp) || isinf(base) || isinf(exp)) {
        fprintf(stderr, "safe_pow: Invalid input (NaN/Inf)\n");
        return 0.0;
    }
    if (base == 0.0 && exp < 0.0) {
        fprintf(stderr, "safe_pow: Division by zero\n");
        return INFINITY;
    }
    errno = 0;
    double result = pow(base, exp);
    if (errno != 0 || isnan(result) || isinf(result)) {
        fprintf(stderr, "safe_pow: Overflow/underflow or error %d\n", errno);
        return (exp > 0.0) ? INFINITY : 0.0;
    }
    return result;
}

/**
 * Safe log: log(base) with checks
 * @param base Base >0
 * @return log(base) or -INFINITY on error
 */
double safe_log(double base) {
    if (base <= 0.0 || isnan(base)) {
        fprintf(stderr, "safe_log: Invalid base <=0 or NaN\n");
        return -INFINITY;
    }
    errno = 0;
    double result = log(base);
    if (errno != 0 || isnan(result)) {
        fprintf(stderr, "safe_log: Error %d\n", errno);
        return -INFINITY;
    }
    return result;
}

/**
 * Allocate 1D double array with error check
 * @param size Size
 * @return Pointer or NULL on failure
 */
double *alloc_darray(size_t size) {
    if (size == 0) {
        fprintf(stderr, "alloc_darray: Zero size\n");
        return NULL;
    }
    double *arr = (double *) malloc(size * sizeof(double));
    if (arr == NULL) {
        fprintf(stderr, "alloc_darray: Malloc failed for %zu elements\n", size);
        return NULL;
    }
    return arr;
}

/**
 * Free array and set pointer to NULL (safe free)
 * @param ptr Pointer
 * @return Void
 */
void free_darray(double **ptr) {
    if (ptr && *ptr) {
        free(*ptr);
        *ptr = NULL;
    }
}

/**
 * Unit conversion: Msun to kg
 * @param mass_Msun Mass [Msun]
 * @return Mass [kg]
 */
double msun_to_kg(double mass_Msun) {
    return mass_Msun * 1.9885e30;
}

/**
 * Unit conversion: kpc to m
 * @param dist_kpc Distance [kpc]
 * @return Distance [m]
 */
double kpc_to_m(double dist_kpc) {
    return dist_kpc * 3.08568e19;
}

/**
 * Unit conversion: km/s to m/s
 * @param vel_kms Velocity [km/s]
 * @return Velocity [m/s]
 */
double kms_to_ms(double vel_kms) {
    return vel_kms * 1000.0;
}

/**
 * Error handling macro wrapper (e.g., for cosmology computations)
 * Log error and exit if fatal
 * @param cond Condition
 * @param msg Error message
 * @return Void (exits on error)
 */
void check_error(int cond, const char *msg) {
    if (cond) {
        fprintf(stderr, "ERROR: %s\n", msg);
        exit(EXIT_FAILURE);
    }
}

/**
 * Cosmology: Comoving distance chi(z) approx integral c/H(z) dz (numerical stub)
 * Use integration.h for full; here simplified for z<<1
 * @param z Redshift
 * @param h h value
 * @return chi(z) [Mpc]
 */
double comoving_distance(double z, double h) {
    double H0_full = 100.0 * h;  /* km/s/Mpc */
    double Hz = H0_full * sqrt(OMEGA_M * pow(1.0 + z, 3.0) + OMEGA_L);  /* Flat LCDM approx */
    return (3e5 / Hz) * log(1.0 + z);  /* c=3e5 km/s; rough for low z */
}

/* Additional: String utils, file I/O helpers if needed */
/**
 * Safe file open
 * @param path Path
 * @param mode Mode ('r','w')
 * @return FILE* or NULL
 */
FILE *safe_fopen(const char *path, const char *mode) {
    FILE *f = fopen(path, mode);
    if (f == NULL) {
        fprintf(stderr, "safe_fopen: Cannot open %s in mode %s: %s\n", path, mode, strerror(errno));
    }
    return f;
}