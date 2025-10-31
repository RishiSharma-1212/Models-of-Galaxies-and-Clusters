/**
 * density.c
 *
 * Core functions for computing density profiles in galaxy and cluster models.
 * Implements NFW (DM halo), Beta (ICM gas), Hernquist (spheroid), Plummer (cored).
 * Computes rho(r), M(<r), Phi(r) with analytic formulas where possible.
 *
 * Units: r in kpc, rho/M in Msun/kpc^3, Phi in (km/s)^2 (G=4.302e-3 (km/s)^2 kpc Msun^{-1}).
 * Spherical symmetry assumed; for arrays, use for loops or vectorize externally.
 *
 * References:
 * - NFW: Navarro et al. (1996)
 * - Beta: Cavaliere & Fusco-Femiano (1976)
 * - Hernquist: Hernquist (1990)
 * - Plummer: Plummer (1911)
 *
 * Usage: Include "density.h"; link with -lm.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "density.h"  // Assumed header with prototypes, constants

/* Constants */
#define PI 3.14159265358979323846
#define G 4.302e-3  /* km^2 s^{-2} kpc Msun^{-1} */

/* Global parameters (editable or passed via struct) */
static double rho_crit = 1.0e9;  // Critical density Msun/kpc^3 (cosmology-dependent)

/**
 * NFW density rho(r) = rho_s / [(r/rs) (1 + r/rs)^2]
 * @param r Radius [kpc]
 * @param rs Scale radius [kpc]
 * @param rho_s Char density [Msun/kpc^3]
 * @return rho(r) [Msun/kpc^3]
 */
double nfw_density(double r, double rs, double rho_s) {
    if (r <= 0.0 || rs <= 0.0) return 0.0;
    double x = r / rs;
    return rho_s / (x * pow(1.0 + x, 2.0));
}

/**
 * NFW enclosed mass M(<r) = 4 pi rho_s rs^3 [ln(1+x) - x/(1+x)], x=r/rs
 * @param r Radius [kpc]
 * @param rs Scale radius [kpc]
 * @param rho_s Char density [Msun/kpc^3]
 * @return M(<r) [Msun]
 */
double nfw_enclosed_mass(double r, double rs, double rho_s) {
    if (r <= 0.0 || rs <= 0.0) return 0.0;
    double x = r / rs;
    double log_term = log(1.0 + x);
    double frac_term = x / (1.0 + x);
    return 4.0 * PI * rho_s * rs * rs * rs * (log_term - frac_term);
}

/**
 * NFW potential Phi(r) = -4 pi G rho_s rs^2 ln(1+x)/x
 * @param r Radius [kpc]
 * @param rs Scale radius [kpc]
 * @param rho_s Char density [Msun/kpc^3]
 * @return Phi(r) [(km/s)^2]
 */
double nfw_potential(double r, double rs, double rho_s) {
    if (r <= 0.0 || rs <= 0.0) return 0.0;
    double x = r / rs;
    return -4.0 * PI * G * rho_s * rs * rs * log(1.0 + x) / x;
}

/**
 * Beta model density rho(r) = rho_0 [1 + (r/rc)^2]^(-3 beta / 2)
 * @param r Radius [kpc]
 * @param rc Core radius [kpc]
 * @param beta Shape param (~0.67)
 * @param rho_0 Central density [Msun/kpc^3]
 * @return rho(r) [Msun/kpc^3]
 */
double beta_density(double r, double rc, double beta, double rho_0) {
    if (r <= 0.0 || rc <= 0.0) return 0.0;
    double x = r / rc;
    return rho_0 * pow(1.0 + x * x, -1.5 * beta);
}

/**
 * Beta enclosed mass: Numerical integration stub (4 pi int r^2 rho dr); analytic for beta !=1/3
 * Simplified: Use quad or approximate; here, return placeholder (implement quad externally)
 * Full analytic involves hypergeometric; for demo, approximate M(<r) ~ 4 pi rho_0 rc^3 * (x^3 / (3 beta - 0.5)) for small x
 * @param r Radius [kpc]
 * @param rc Core radius [kpc]
 * @param beta Shape param
 * @param rho_0 Central density [Msun/kpc^3]
 * @return M(<r) [Msun] (approximate)
 */
double beta_enclosed_mass(double r, double rc, double beta, double rho_0) {
    if (r <= 0.0 || rc <= 0.0) return 0.0;
    double x = r / rc;
    // Approximate for beta=2/3: M(<r) = 4 pi rho_0 rc^3 * x^3 / (1 + x^2)^{3 beta / 2} * factor
    // Better: Use integral; here simplified
    return 4.0 * PI * rho_0 * rc * rc * rc * (x * x * x) / (3.0 * beta - 0.5);  // Placeholder; use quad for accuracy
}

/**
 * Hernquist density rho(r) = (M a) / (2 pi r (r + a)^3)
 * @param r Radius [kpc]
 * @param a Scale [kpc]
 * @param M Total mass [Msun]
 * @return rho(r) [Msun/kpc^3]
 */
double hernquist_density(double r, double a, double M) {
    if (r <= 0.0 || a <= 0.0) return 0.0;
    return (M * a) / (2.0 * PI * r * pow(r + a, 3.0));
}

/**
 * Hernquist enclosed mass M(<r) = M (r / (r + a))^2
 * @param r Radius [kpc]
 * @param a Scale [kpc]
 * @param M Total mass [Msun]
 * @return M(<r) [Msun]
 */
double hernquist_enclosed_mass(double r, double a, double M) {
    if (r <= 0.0 || a <= 0.0) return 0.0;
    return M * pow(r / (r + a), 2.0);
}

/**
 * Plummer density rho(r) = (3 M / 4 pi a^3) (1 + (r/a)^2)^{-5/2}
 * @param r Radius [kpc]
 * @param a Scale [kpc]
 * @param M Total mass [Msun]
 * @return rho(r) [Msun/kpc^3]
 */
double plummer_density(double r, double a, double M) {
    if (r <= 0.0 || a <= 0.0) return 0.0;
    double x = r / a;
    double pref = 3.0 * M / (4.0 * PI * pow(a, 3.0));
    return pref * pow(1.0 + x * x, -2.5);
}

/**
 * Plummer enclosed mass M(<r) = M r^3 / (a^2 + r^2)^{3/2}
 * @param r Radius [kpc]
 * @param a Scale [kpc]
 * @param M Total mass [Msun]
 * @return M(<r) [Msun]
 */
double plummer_enclosed_mass(double r, double a, double M) {
    if (r <= 0.0 || a <= 0.0) return 0.0;
    double denom = pow(a * a + r * r, 1.5);
    return M * pow(r, 3.0) / denom;
}

/* Vectorized versions for arrays (N points) */
void nfw_density_vec(double *r, int N, double rs, double rho_s, double *rho_out) {
    for (int i = 0; i < N; i++) {
        rho_out[i] = nfw_density(r[i], rs, rho_s);
    }
}

void nfw_enclosed_mass_vec(double *r, int N, double rs, double rho_s, double *m_out) {
    for (int i = 0; i < N; i++) {
        m_out[i] = nfw_enclosed_mass(r[i], rs, rho_s);
    }
}

/* Similar vec functions for others... (omitted for brevity; implement analogously) */

/* Normalization helpers, e.g., rho_s for NFW given M_vir, c */
double compute_nfw_rho_s(double M_vir, double c, double rho_crit) {
    double f_c = log(1.0 + c) - c / (1.0 + c);
    return (M_vir * (200.0 / 3.0)) / (4.0 * PI * pow((M_vir / (200.0 * rho_crit)), 1.0/3.0) * f_c);
    // rs = R_vir / c, R_vir = (3 M_vir / (4 pi 200 rho_crit))^{1/3}
}