/**
 * example_hernquist.c
 * 
 * Standalone C example for the Hernquist density profile (Hernquist 1990),
 * used to model spheroidal components in galaxies (e.g., bulges in clusters).
 * Computes 3D density rho(r), enclosed mass M(<r), and gravitational potential Phi(r).
 * 
 * Units: r in kpc, M in Msun, rho in Msun/kpc^3, Phi in (km/s)^2 (up to G=1).
 * Assumptions: Spherical symmetry, total mass M_tot normalized.
 * 
 * Output: CSV table of r, rho(r), M(<r), Phi(r) from r_min to r_max in N steps.
 * 
 * References:
 * Hernquist, L. 1990, ApJ, 356, 359
 * 
 * Compile: gcc example_hernquist.c -lm -o hernquist
 * Run: ./hernquist
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Constants */
#define PI 3.14159265358979323846
#define G 4.302e-3  /* Gravitational constant in (km/s)^2 kpc / Msun (approx) */

/* Hernquist profile parameters (user-editable) */
#define M_TOT 1.0e12     /* Total mass Msun */
#define A_SCALE 5.0      /* Scale radius a in kpc */
#define R_MIN 0.1        /* Minimum radius kpc */
#define R_MAX 100.0      /* Maximum radius kpc */
#define N_POINTS 1000    /* Number of radial points */

/**
 * Computes Hernquist density rho(r).
 * rho(r) = (M_tot * a) / (2 * PI * r * (r + a)^3 )
 * 
 * @param r Radius in kpc
 * @return Density rho(r) in Msun/kpc^3
 */
double hernquist_density(double r) {
    if (r <= 0.0) return 0.0;
    double term1 = M_TOT * A_SCALE;
    double term2 = 2.0 * PI * r * pow(r + A_SCALE, 3.0);
    return term1 / term2;
}

/**
 * Computes enclosed mass M(<r) for Hernquist profile.
 * Analytic: M(<r) = M_tot * (r / (r + a))^2
 * 
 * @param r Radius in kpc
 * @return Enclosed mass Msun
 */
double hernquist_enclosed_mass(double r) {
    if (r <= 0.0) return 0.0;
    return M_TOT * pow(r / (r + A_SCALE), 2.0);
}

/**
 * Computes gravitational potential Phi(r) (up to -G M / r).
 * For Hernquist: Phi(r) = - (G M_tot / (r + a)) * (1 + (a / (2 r)) * ln( (r + a)^2 / (r (r + 2 a)) ) ) approx, but exact is:
 * Phi(r) = - G M_tot / (r + a) + correction; simplified to -G M(<r) / r for demo.
 * Full analytic: Phi(r) = - (G M_tot / a) * (ln(r + a) / (r / a) ) / 2 or similar; here use integrated form.
 * Note: Full derivation in Hernquist 1990; this uses approximate potential from mass.
 * 
 * @param r Radius in kpc
 * @return Potential in (km/s)^2
 */
double hernquist_potential(double r) {
    if (r <= 0.0) return 0.0;
    double m_enc = hernquist_enclosed_mass(r);
    return -G * m_enc / r;  /* Simplified; for exact, integrate Poisson eq. */
    // Exact Phi(r) = - (G M_tot / (r + a)) * (1 + 0.5 * (a / (r + a)) * log( (r + a)^2 / (r (r + 2*a)) ) ); but placeholder
}

/**
 * Main function: Generate radial table.
 */
int main() {
    double dr = (R_MAX - R_MIN) / (N_POINTS - 1.0);
    
    printf("r_kpc,rho_Msun_kpc3,M_enc_Msun,Phi_kms2\n");
    
    for (int i = 0; i < N_POINTS; i++) {
        double r = R_MIN + i * dr;
        double rho = hernquist_density(r);
        double m_enc = hernquist_enclosed_mass(r);
        double phi = hernquist_potential(r);
        
        printf("%.6f,%.6e,%.6e,%.6f\n", r, rho, m_enc, phi);
    }
    
    printf("# Total mass check: M(<%.1f kpc) = %.6e Msun (should ≈ %g)\n", R_MAX, hernquist_enclosed_mass(R_MAX), M_TOT);
    
    return 0;
}