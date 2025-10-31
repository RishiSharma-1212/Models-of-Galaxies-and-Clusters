/**
 * density.h
 *
 * Header file for density.c: Density profiles for galaxy and cluster models.
 * Declares functions for NFW (DM), Beta (gas), Hernquist (spheroid), Plummer (cored).
 * Includes scalar and vectorized (array) versions for efficiency.
 *
 * Units: r [kpc], rho/M [Msun/kpc^3], Phi [(km/s)^2] with G=4.302e-3.
 *
 * Usage: #include "density.h"; compile with -lm.
 *
 * References:
 * - NFW: Navarro et al. (1996) [web:159][web:165]
 * - Beta: Cavaliere & Fusco-Femiano (1976)
 * - Hernquist: Hernquist (1990) [web:189][web:194]
 * - Plummer: Plummer (1911) [web:212]
 */

#ifndef DENSITY_H
#define DENSITY_H

/* Constants */
#define PI 3.14159265358979323846
#define G 4.302e-3  /* km^2 s^{-2} kpc Msun^{-1} (astrophysical units) */
#define RHO_CRIT_DEFAULT 1.878e-29  /* Default critical density [h^2 Msun / kpc^3] for h=0.7 */

/* Struct for profile parameters (optional, for passing configs) */
typedef struct {
    double mass;     /* Total/virial mass [Msun] */
    double scale_r;  /* Scale radius (rs, rc, a) [kpc] */
    double rho_char; /* Characteristic density [Msun/kpc^3] */
    double beta;     /* For Beta model */
    double conc;     /* Concentration c for NFW */
} ProfileParams;

/**
 * NFW Functions
 */

/**
 * NFW density: rho(r) = rho_s / [(r/rs) (1 + r/rs)^2]
 * @param r Radius [kpc]
 * @param rs Scale radius [kpc]
 * @param rho_s Characteristic density [Msun/kpc^3]
 * @return rho(r) [Msun/kpc^3]
 */
double nfw_density(double r, double rs, double rho_s);

/**
 * NFW enclosed mass: M(<r) = 4 pi rho_s rs^3 [ln(1+x) - x/(1+x)], x=r/rs
 * @param r Radius [kpc]
 * @param rs Scale radius [kpc]
 * @param rho_s Characteristic density [Msun/kpc^3]
 * @return M(<r) [Msun]
 */
double nfw_enclosed_mass(double r, double rs, double rho_s);

/**
 * NFW potential: Phi(r) = -4 pi G rho_s rs^2 ln(1+x)/x
 * @param r Radius [kpc]
 * @param rs Scale radius [kpc]
 * @param rho_s Characteristic density [Msun/kpc^3]
 * @return Phi(r) [(km/s)^2]
 */
double nfw_potential(double r, double rs, double rho_s);

/**
 * Vectorized NFW density for array r[0..N-1]
 * @param r Input array of radii [kpc]
 * @param N Number of points
 * @param rs Scale radius [kpc]
 * @param rho_s Characteristic density [Msun/kpc^3]
 * @param rho_out Output array for densities
 */
void nfw_density_vec(double *r, int N, double rs, double rho_s, double *rho_out);

/**
 * Vectorized NFW enclosed mass
 */
void nfw_enclosed_mass_vec(double *r, int N, double rs, double rho_s, double *m_out);

/**
 * Beta Model Functions
 */

/**
 * Beta density: rho(r) = rho_0 [1 + (r/rc)^2]^(-3 beta / 2)
 * @param r Radius [kpc]
 * @param rc Core radius [kpc]
 * @param beta Shape parameter (~0.67)
 * @param rho_0 Central density [Msun/kpc^3]
 * @return rho(r) [Msun/kpc^3]
 */
double beta_density(double r, double rc, double beta, double rho_0);

/**
 * Beta enclosed mass (approximate; use numerical integration for precision)
 * @param r Radius [kpc]
 * @param rc Core radius [kpc]
 * @param beta Shape parameter
 * @param rho_0 Central density [Msun/kpc^3]
 * @return M(<r) [Msun]
 */
double beta_enclosed_mass(double r, double rc, double beta, double rho_0);

/**
 * Hernquist Functions
 */

/**
 * Hernquist density: rho(r) = (M a) / (2 pi r (r + a)^3)
 * @param r Radius [kpc]
 * @param a Scale radius [kpc]
 * @param M Total mass [Msun]
 * @return rho(r) [Msun/kpc^3]
 */
double hernquist_density(double r, double a, double M);

/**
 * Hernquist enclosed mass: M(<r) = M (r / (r + a))^2
 * @param r Radius [kpc]
 * @param a Scale radius [kpc]
 * @param M Total mass [Msun]
 * @return M(<r) [Msun]
 */
double hernquist_enclosed_mass(double r, double a, double M);

/**
 * Plummer Functions
 */

/**
 * Plummer density: rho(r) = (3 M / 4 pi a^3) (1 + (r/a)^2)^{-5/2}
 * @param r Radius [kpc]
 * @param a Scale radius [kpc]
 * @param M Total mass [Msun]
 * @return rho(r) [Msun/kpc^3]
 */
double plummer_density(double r, double a, double M);

/**
 * Plummer enclosed mass: M(<r) = M r^3 / (a^2 + r^2)^{3/2}
 * @param r Radius [kpc]
 * @param a Scale radius [kpc]
 * @param M Total mass [Msun]
 * @return M(<r) [Msun]
 */
double plummer_enclosed_mass(double r, double a, double M);

/**
 * Normalization helper: Compute rho_s for NFW given M_vir, c, rho_crit
 * @param M_vir Virial mass [Msun]
 * @param c Concentration
 * @param rho_crit Critical density [Msun/kpc^3]
 * @return rho_s [Msun/kpc^3]
 */
double compute_nfw_rho_s(double M_vir, double c, double rho_crit);

#endif /* DENSITY_H */