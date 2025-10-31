/**
 * profiles.c
 *
 * Higher-level functions for full galaxy/cluster profiles, combining DM (NFW),
 * gas (Beta), and stellar (Hernquist/Plummer) components.
 * Computes total density rho_tot(r), total enclosed mass M_tot(<r) (analytic + num),
 * and pressure/temperature for hydrostatic equilibrium.
 *
 * Uses density.h and integration.h; assumes spherical symmetry.
 * Sample: Cluster with M_vir=1e15 Msun, f_gas=0.15, f_stars=0.01.
 *
 * Units: r [kpc], rho/M [Msun/kpc^3], P [keV cm^{-3}] for gas.
 *
 * References:
 * - Combined profiles: NFW+Beta for clusters [web:153][web:326][web:327]
 * - Gas/stellar fractions: Observations [web:157][web:334]
 *
 * Usage: #include "profiles.h"; link with density.o integration.o -lm.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "density.h"
#include "integration.h"
#include "profiles.h"  // Assumed header with prototypes

/* Constants */
#define PI 3.14159265358979323846
#define K_B 8.617e-5   /* Boltzmann keV/K */
#define MU 0.59        /* Mean mol. wt. for ICM */
#define M_P 1.67e-24   /* Proton mass g; use Msun/kpc conversions externally */

/* Sample cluster parameters (editable via struct) */
static double M_VIR = 1.0e15;      /* Virial mass Msun */
static double CONC = 8.0;          /* DM concentration */
static double F_GAS = 0.15;        /* Gas fraction */
static double F_STARS = 0.01;      /* Stellar fraction */
static double T_GAS = 5.0e7;       /* Gas temperature K (~5 keV) */
static double RHO_CRIT = 2.775e11; /* h^2 Msun kpc^{-3} approx */

/* Cluster config struct */
typedef struct {
    double m_vir;     /* Virial mass Msun */
    double conc;      /* NFW c */
    double f_gas;     /* Gas frac */
    double f_stars;   /* Stars frac */
    double t_gas;     /* Gas T [K] */
    double rho_crit;  /* Critical density */
} ClusterParams;

/**
 * Initialize cluster params (default or custom).
 * @param params Output struct
 * @param m_vir Optional M_vir
 * @param conc Optional c
 * etc.
 */
void init_cluster_params(ClusterParams *params, double m_vir, double conc, double f_gas, double f_stars, double t_gas, double rho_crit) {
    params->m_vir = (m_vir > 0.0) ? m_vir : M_VIR;
    params->conc = (conc > 0.0) ? conc : CONC;
    params->f_gas = (f_gas > 0.0) ? f_gas : F_GAS;
    params->f_stars = (f_stars > 0.0) ? f_stars : F_STARS;
    params->t_gas = (t_gas > 0.0) ? t_gas : T_GAS;
    params->rho_crit = (rho_crit > 0.0) ? rho_crit : RHO_CRIT;
}

/**
 * Total density rho_tot(r) = rho_DM(r) + rho_gas(r) + rho_stars(r)
 * NFW DM, Beta gas, Hernquist stars.
 * @param r Radius [kpc]
 * @param params Cluster params
 * @return rho_tot(r) [Msun/kpc^3]
 */
double cluster_total_density(double r, ClusterParams *params) {
    if (r <= 0.0) return 0.0;

    /* DM: NFW */
    double rs_dm = (3.0 * params->m_vir / (4.0 * PI * 200.0 * params->rho_crit * pow(params->conc, 3.0))) * params->conc;  /* Simplified rs */
    double rho_s_dm = compute_nfw_rho_s(params->m_vir, params->conc, params->rho_crit);
    double rho_dm = nfw_density(r, rs_dm, rho_s_dm);

    /* Gas: Beta, normalized to f_gas M_vir */
    double rc_gas = 0.1 * rs_dm;  /* Core ~10% rs */
    double beta_gas = 2.0 / 3.0;
    double m_gas = params->f_gas * params->m_vir;
    double rho_0_gas = m_gas / (4.0 * PI * pow(rc_gas, 3.0) * (1.0 / (3.0 * beta_gas - 0.5)));  /* Approx normalization */
    double rho_gas = beta_density(r, rc_gas, beta_gas, rho_0_gas);

    /* Stars: Hernquist for BCG, f_stars M_vir */
    double a_stars = 0.01 * rs_dm;  /* Compact bulge */
    double m_stars = params->f_stars * params->m_vir;
    double rho_stars = hernquist_density(r, a_stars, m_stars);

    return rho_dm + rho_gas + rho_stars;
}

/**
 * Total enclosed mass M_tot(<r) = M_DM(<r) + M_gas(<r) + M_stars(<r)
 * Uses analytic where possible, numerical fallback.
 * @param r Radius [kpc]
 * @param params Cluster params
 * @param use_numeric Flag for numerical integration
 * @return M_tot(<r) [Msun]
 */
double cluster_total_mass(double r, ClusterParams *params, int use_numeric) {
    if (r <= 0.0) return 0.0;

    /* DM: NFW analytic or num */
    double rs_dm = ...;  /* As above */
    double rho_s_dm = ...;
    double m_dm = use_numeric ? numerical_nfw_mass(r, rs_dm, rho_s_dm, 1000) : nfw_enclosed_mass(r, rs_dm, rho_s_dm);
