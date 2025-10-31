/**
 * profiles.h
 *
 * Header file for profiles.c: Combined galaxy/cluster profiles (DM + gas + stars).
 * Declares structs and prototypes for total density/mass, gas pressure,
 * integrating NFW (DM), Beta (gas), Hernquist/Plummer (stars).
 *
 * Assumes spherical symmetry; uses density.h and integration.h for components.
 * Sample: Cluster with f_gas=0.15, f_stars=0.01, isothermal gas.
 *
 * Units: r [kpc], rho/M [Msun/kpc^3], P [keV cm^{-3}] approx.
 *
 * References:
 * - NFW+Beta clusters: Lokas (2001); Voit et al. (2005) [web:153][web:341]
 * - Combined models: Observations and simulations [web:327][web:334][web:342]
 *
 * Usage: #include "profiles.h"; #include "density.h"; #include "integration.h"; link -lm.
 */

#ifndef PROFILES_H
#define PROFILES_H

#include "density.h"  /* For base profiles */

/* Constants */
#define PI 3.14159265358979323846
#define K_B 8.617e-5   /* Boltzmann constant keV/K */
#define MU_GAS 0.59    /* Mean molecular weight for ICM (ionized H/He) */
#define M_PROTON 1.67e-24  /* g; convert to Msun/kpc^3 as needed */
#define RHO_CRIT_DEFAULT 2.775e11  /* h^2 Msun kpc^{-3} approx (h=0.7) */

/* Cluster configuration struct */
typedef struct {
    double m_vir;      /* Virial mass Msun */
    double conc;       /* NFW concentration c */
    double f_gas;      /* Gas mass fraction (M_gas / M_vir) */
    double f_stars;    /* Stellar mass fraction (M_stars / M_vir) */
    double t_gas;      /* Gas temperature [K] (isothermal) */
    double rho_crit;   /* Critical density [Msun kpc^{-3}] */
    double rc_gas;     /* Gas core radius [kpc] (relative to rs_dm) */
    double a_stars;    /* Stellar scale radius [kpc] (Hernquist) */
    double beta_gas;   /* Beta model slope (~2/3) */
} ClusterParams;

/**
 * Initialize cluster parameters (defaults or custom).
 * @param params Output: ClusterParams struct
 * @param m_vir Virial mass [Msun] (default 1e15)
 * @param conc NFW concentration (default 8.0)
 * @param f_gas Gas fraction (default 0.15)
 * @param f_stars Stellar fraction (default 0.01)
 * @param t_gas Gas temperature [K] (default 5e7 ~5 keV)
 * @param rho_crit Critical density [Msun kpc^{-3}] (default 2.775e11)
 * @return Void
 */
void init_cluster_params(ClusterParams *params, double m_vir, double conc, 
                         double f_gas, double f_stars, double t_gas, double rho_crit);

/**
 * Total density rho_tot(r) = rho_NFW(r) + rho_Beta(r) + rho_Hernquist(r)
 * @param r Radius [kpc]
 * @param params Cluster parameters
 * @return rho_tot(r) [Msun kpc^{-3}]
 */
double cluster_total_density(double r, ClusterParams *params);

/**
 * Total enclosed mass M_tot(<r) = M_NFW(<r) + M_Beta(<r) + M_Hernquist(<r)
 * Analytic where possible, numerical fallback.
 * @param r Radius [kpc]
 * @param params Cluster parameters
 * @param use_numeric Flag: 1 for numerical integration, 0 for analytic
 * @return M_tot(<r) [Msun]
 */
double cluster_total_mass(double r, ClusterParams *params, int use_numeric);

/**
 * Gas pressure P(r) = [rho_gas(r) / (mu m_p)] kT (isothermal)
 * @param r Radius [kpc]
 * @param params Cluster parameters
 * @return P(r) [keV cm^{-3}] (approximate; unit conversion assumed)
 */
double cluster_gas_pressure(double r, ClusterParams *params);

/**
 * Vectorized total density: rho_tot(r[i]) for array
 * @param r Input array of radii [kpc]
 * @param N Number of points
 * @param params Cluster parameters
 * @param rho_out Output array of densities [Msun kpc^{-3}]
 */
void cluster_total_density_vec(double *r, int N, ClusterParams *params, double *rho_out);

/**
 * Verify total mass at R_vir ≈ M_vir (relative error)
 * @param params Cluster parameters
 * @return Relative error (fraction)
 */
double verify_total_mass(ClusterParams *params);

/**
 * Projected surface density Sigma(R) = 2 ∫ rho_tot(sqrt(R^2 + l^2)) dl (l.o.s.)
 * Stub: Requires numerical integration over line-of-sight; placeholder.
 * @param R Projected radius [kpc]
 * @param params Cluster parameters
 * @return Sigma(R) [Msun kpc^{-2}]
 */
double projected_density(double R, ClusterParams *params);

/* Optional extensions: e.g., for Plummer stars or temperature profile */
/* double cluster_stellar_density(double r, ClusterParams *params); */

#endif /* PROFILES_H */
