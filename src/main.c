/**
 * main.c
 *
 * Entry point for Models-of-Galaxies-and-Clusters: Demonstrates density profiles
 * and numerical integration for a sample galaxy cluster halo.
 * Computes NFW density rho(r), enclosed mass M(<r) (analytic + numerical), and potential.
 * Outputs radial table to stdout (CSV format) for plotting (e.g., with Gnuplot).
 *
 * Sample: M_vir=1e15 Msun cluster at z=0, c=10, using Planck cosmology (rho_crit approx).
 *
 * Usage: gcc main.c density.o integration.o -lm -o main; ./main > output.csv
 *
 * References:
 * - NFW: Navarro et al. (1996) [web:165][web:314]
 * - Simulations: Initial conditions and profiles [web:315][web:317]
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "density.h"
#include "integration.h"

/* Constants for sample cosmology (Planck18 approx, z=0) */
#define RHO_CRIT 2.775e11  /* h^2 Msun h^{-1} Mpc^{-3} -> Msun/kpc^3 (h=0.674) */
#define H 0.674  /* Hubble const */
#define OMEGA_M 0.315  /* Matter density */

/* Sample cluster parameters */
#define M_VIR 1.0e15      /* Virial mass Msun */
#define CONC 10.0         /* Concentration c */
#define Z_INIT 0.0        /* Initial redshift (z=0 for simplicity) */
#define R_MIN 0.1         /* Min radius kpc */
#define R_MAX 1500.0      /* Approx R_vir ~ (3 M_vir / (4 pi 200 rho_crit))^{1/3} */
#define N_POINTS 1000     /* Radial points */

/**
 * Main function: Compute and output NFW profile table.
 */
int main() {
    /* Compute NFW parameters */
    double rho_crit = RHO_CRIT * pow(H, 2.0);  /* Adjusted for units */
    double rs = (3.0 * M_VIR / (4.0 * PI * 200.0 * rho_crit)) / CONC;  /* Scale radius kpc */
    double rho_s = compute_nfw_rho_s(M_VIR, CONC, rho_crit);  /* Char density */
    double r_vir = rs * CONC;  /* Virial radius kpc */

    /* Radial grid */
    double dr = (R_MAX - R_MIN) / (N_POINTS - 1.0);
    double *r = (double*) malloc(N_POINTS * sizeof(double));
    for (int i = 0; i < N_POINTS; i++) {
        r[i] = R_MIN + i * dr;
    }

    /* Output header */
    printf("r_kpc,rho_Msun_kpc3,M_enc_analytic_Msun,M_enc_numerical_Msun,Phi_kms2\n");

    /* Compute and print for each r */
    for (int i = 0; i < N_POINTS; i++) {
        double rr = r[i];
        double rho = nfw_density(rr, rs, rho_s);
        double m_analytic = nfw_enclosed_mass(rr, rs, rho_s);
        double m_numerical = numerical_nfw_mass(rr, rs, rho_s, 1000);  /* N=1000 subintervals */
        double phi = nfw_potential(rr, rs, rho_s);

        /* CSV output */
        printf("%.6f,%.6e,%.6e,%.6e,%.6f\n", rr, rho, m_analytic, m_numerical, phi);

        /* Progress every 100 points */
        if ((i + 1) % 100 == 0) {
            fprintf(stderr, "Processed %d/%d points (r=%.1f kpc)\n", i + 1, N_POINTS, rr);
        }
    }

    /* Verification at R_vir */
    double m_vir_check = nfw_enclosed_mass(r_vir, rs, rho_s);
    printf("# Verification: M(<%.1f kpc) = %.6e Msun (target: %g Msun)\n", 
           r_vir, m_vir_check, M_VIR);
    double rel_error = fabs((m_vir_check - M_VIR) / M_VIR) * 100.0;
    printf("# Relative error: %.4f%%\n", rel_error);

    /* Free memory */
    free(r);

    return 0;
}

/* Optional: Extend for other profiles, e.g., Beta gas or Hernquist BCG */
/* Example stub:
   // For Beta: double rc=100.0, beta=0.67, rho_0=1e10;
   // rho = beta_density(rr, rc, beta, rho_0);
   // m_num = numerical_enclosed_mass(rr, 1000, beta_rho_wrapper, &beta_params);
*/