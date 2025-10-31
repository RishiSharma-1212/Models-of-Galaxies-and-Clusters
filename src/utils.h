/**
 * utils.h
 *
 * Header file for utils.c: Astronomy utilities for cosmology, safe math,
 * array management, unit conversions, and error handling.
 * Declares functions for rho_crit, safe pow/log, alloc/free, and basic cosmology.
 *
 * Cosmology: Defaults to Planck18 (h=0.674, Om=0.315); functions take optional params.
 * Safe math: Checks for NaN/Inf/overflow in pow/log.
 * Arrays: Dynamic double arrays with error handling.
 *
 * Units: Focus on astrophysical (kpc, Msun, km/s); conversions to SI.
 *
 * References:
 * - Cosmology: Planck Collaboration (2020); Hogg (1999) [web:388][web:389]
 * - Math: Numerical Recipes (Press et al. 2007) [web:374]
 * - Constants: Astropy/yt [web:369][web:387][web:393]
 *
 * Usage: #include "utils.h"; compile with -lm.
 */

#ifndef UTILS_H
#define UTILS_H

#include <stddef.h>  /* For size_t */

/* Cosmology constants (Planck18 defaults) */
#define H0_DEFAULT 67.4          /* km/s/Mpc */
#define OMEGA_M_DEFAULT 0.315    /* Matter density */
#define OMEGA_L_DEFAULT 0.685    /* Lambda */
#define OMEGA_B_DEFAULT 0.049    /* Baryon */
#define H_DEFAULT 0.674          /* h = H0/100 */

/* Physical constants */
#define MSUN_KG 1.9885e30        /* Solar mass [kg] */
#define KPC_M 3.08568e19         /* kpc [m] */
#define G_SI 6.67430e-11         /* Grav const m^3 kg^{-1} s^{-2} */

/**
 * Compute critical density rho_crit(z=0) = 3 H^2 / (8 pi G) [Msun kpc^{-3}]
 * @param h Optional h (default H_DEFAULT)
 * @return rho_crit [Msun kpc^{-3}]
 */
double compute_rho_crit(double h);

/**
 * Safe power: pow(base, exp) with NaN/Inf/overflow checks
 * @param base Base value
 * @param exp Exponent
 * @return pow(base, exp) or 0.0/Inf on error
 */
double safe_pow(double base, double exp);

/**
 * Safe natural log: log(base) with base >0 check
 * @param base Base >0
 * @return log(base) or -Inf on error
 */
double safe_log(double base);

/**
 * Allocate 1D double array
 * @param size Number of elements
 * @return Pointer to array or NULL on failure
 */
double *alloc_darray(size_t size);

/**
 * Safe free for double array (sets ptr to NULL)
 * @param ptr Pointer to pointer
 * @return Void
 */
void free_darray(double **ptr);

/**
 * Unit conversion: Msun to kg
 * @param mass_Msun [Msun]
 * @return mass [kg]
 */
double msun_to_kg(double mass_Msun);

/**
 * Unit conversion: kpc to m
 * @param dist_kpc [kpc]
 * @return dist [m]
 */
double kpc_to_m(double dist_kpc);

/**
 * Unit conversion: km/s to m/s
 * @param vel_kms [km/s]
 * @return vel [m/s]
 */
double kms_to_ms(double vel_kms);

/**
 * Error check: Exit on failure
 * @param cond Condition (non-zero error)
 * @param msg Error message
 * @return Void (exits if cond true)
 */
void check_error(int cond, const char *msg);

/**
 * Approximate comoving distance chi(z) [Mpc] (low-z; use integration for full)
 * @param z Redshift
 * @param h h parameter
 * @return chi(z) [Mpc]
 */
double comoving_distance(double z, double h);

/**
 * Safe file open with error message
 * @param path File path
 * @param mode Open mode ('r', 'w', etc.)
 * @return FILE* or NULL on failure
 */
FILE *safe_fopen(const char *path, const char *mode);

/* Optional: Set cosmology params (globals for simplicity) */
/**
 * Set Hubble h
 * @param h New h
 * @return Void
 */
void set_h(double h);

/**
 * Set Omega_m
 * @param om Omega_m
 * @return Void
 */
void set_omega_m(double om);

/* Macro for min/max (if not using <math.h>) */
#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

#endif /* UTILS_H */
