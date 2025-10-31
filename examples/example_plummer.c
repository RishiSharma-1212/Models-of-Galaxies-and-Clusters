/**
 * example_plummer.c
 * 
 * Standalone C example for the Plummer density profile (Plummer 1911),
 * used to model globular clusters, dwarf galaxies, or cored stellar systems in simulations.
 * Computes 3D density rho(r), enclosed mass M(<r), and gravitational potential Phi(r).
 * 
 * Units: r in kpc, M in Msun, rho in Msun/kpc^3, Phi in (km/s)^2 (up to G=1).
 * Assumptions: Spherical symmetry, total mass M_tot normalized.
 * 
 * Output: CSV table of r, rho(r), M(<r), Phi(r) from r_min to r_max in N steps.
 * 
 * References:
 * Plummer, H. C. 1911, MNRAS, 71, 460
 * 
 * Compile: gcc example_plummer.c -lm -o plummer
 * Run: ./plummer > output.csv
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Constants */
#define PI 3.14159265358979323846
#define G 4.302e-3  /* Gravitational constant in (km/s)^2 kpc / Msun (approx) */

/* Plummer profile parameters (user-editable) */
#define M_TOT 1.0e10      /* Total mass M_0 in Msun */
#define A_SCALE 1.0       /* Plummer scale radius a in kpc (core radius) */
#define R_MIN 0.01        /* Minimum radius kpc (avoid singularity at 0) */
#define R_MAX 20.0        /* Maximum radius kpc */
#define N_POINTS 1000     /* Number of radial points */

/**
 * Computes Plummer density rho(r).
 * rho(r) = (3 M_0 / (4 π a^3)) * (1 + (r/a)^2)^{-5/2}
 * 
 * @param r Radius in kpc
 * @return Density rho(r) in Msun/kpc^3
 */
double plummer_density(double r) {
    if (r <= 0.0) return 0.0;
    double x = r / A_SCALE;
    double prefactor = 3.0 * M_TOT / (4.0 * PI * pow(A_SCALE, 3.0));
    return prefactor * pow(1.0 + x * x, -2.5);
}

/**
 * Computes enclosed mass M(<r) for Plummer profile.
 * Analytic: M(<r) = M_0 * (r^3 / (a^2 + r^2)^{3/2})
 * 
 * @param r Radius in kpc
 * @return Enclosed mass in Msun
 */
double plummer_enclosed_mass(double r) {
    if (r <= 0.0) return 0.0;
    double x = r / A_SCALE;
    double denom = pow(A_SCALE * A_SCALE + r * r, 1.5);
    return M_TOT * (r * r * r) / denom;
}

/**
 * Computes gravitational potential Phi(r).
 * Analytic: Phi(r) = - G M_0 / sqrt(r^2 + a^2)
 * 
 * @param r Radius in kpc
 * @return Potential in (km/s)^2
 */
double plummer_potential(double r) {
    if (r <= 0.0) return 0.0;
    double rsq = r * r + A_SCALE * A_SCALE;
    return -G * M_TOT / sqrt(rsq);
}

/**
 * Main function: Generate radial table.
 */
int main() {
    double dr = (R_MAX - R_MIN) / (N_POINTS - 1.0);
    
    // CSV header
    printf("r_kpc,rho_Msun_kpc3,M_enc_Msun,Phi_kms2\n");
    
    for (int i = 0; i < N_POINTS; i++) {
        double r = R_MIN + i * dr;
        double rho = plummer_density(r);
        double m_enc = plummer_enclosed_mass(r);
        double phi = plummer_potential(r);
        
        // Format for CSV
        printf("%.6f,%.6e,%.6e,%.6f\n", r, rho, m_enc, phi);
    }
    
    // Verification
    double m_total_check = plummer_enclosed_mass(R_MAX);
    printf("# Total mass check: M(<%.1f kpc) = %.6e Msun (should ≈ %g Msun)\n", 
           R_MAX, m_total_check, M_TOT);
    
    return 0;
}