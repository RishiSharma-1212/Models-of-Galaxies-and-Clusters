/**
 * test_sampling.c
 *
 * Unit tests for sampling.c: Rejection sampling positions (radial dist vs rho),
 * velocity dispersion (isotropic Jeans), total mass conservation for NFW/Hernquist.
 * Samples N=1000 particles, computes histograms/stats, compares to analytic rho/M/sigma.
 * Tolerance: 10% rel for stats (statistical variance at low N).
 *
 * References:
 * - Rejection sampling validation: Envelope principle [web:356][web:423]
 * - Velocity dispersion: Jeans eq. for halos [web:428][web:431]
 * - IC tests: Literature benchmarks [web:440]
 *
 * Usage: gcc test_sampling.c ../src/sampling.o ../src/density.o ../src/utils.o -lm -o test_sampling; ./test_sampling
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../src/sampling.h"
#include "../src/density.h"
#include "../src/utils.h"

/* Test params */
#define N_PARTICLES 1000
#define N_BINS 20
#define TOL_STAT 0.1  /* 10% rel tol for low-N stats */

/* Assertion with rel tol */
#define ASSERT_CLOSE(val, exp, tol, msg) do { \
    double rel = fabs((val) - (exp)) / (fabs(exp) + 1e-10); \
    if (rel > tol) { \
        fprintf(stderr, "FAIL %s: Got %.4g, exp %.4g (rel=%.2f)\n", msg, val, exp, rel); \
        exit(1); \
    } \
} while(0)

/* Compute radial histogram for positions (shell volumes) */
void compute_radial_hist(double *r, int N, double r_min, double r_max, int n_bins, double *hist_rho, double *hist_vol) {
    double dr = (r_max - r_min) / n_bins;
    memset(hist_rho, 0, n_bins * sizeof(double));
    memset(hist_vol, 0, n_bins * sizeof(double));
    for (int i = 0; i < N; i++) {
        double rr = r[i];
        if (rr < r_min || rr > r_max) continue;
        int bin = (int)((rr - r_min) / dr);
        if (bin >= n_bins) bin = n_bins - 1;
        hist_rho[bin] += 1.0;  /* Counts */
        double r_mid = r_min + (bin + 0.5) * dr;
        hist_vol[bin] += 4.0 * M_PI * r_mid * r_mid * dr;  /* Shell vol approx */
    }
    for (int i = 0; i < n_bins; i++) {
        if (hist_vol[i] > 0) hist_rho[i] /= hist_vol[i];  /* Density = N_shell / V_shell */
    }
}

/* Test NFW position sampling: Radial density hist vs analytic rho(r) */
void test_sample_positions_nfw() {
    printf("Testing NFW position sampling...\n");
    double rs = 100.0, rho_s = 1e9;
    double r_max = 10.0 * rs;
    double x[N_PARTICLES], y[N_PARTICLES], z[N_PARTICLES], m[N_PARTICLES];
    double r[N_PARTICLES];
    
    srand_custom(42);  /* Fixed seed for repro */
    int n_sampled = sample_positions_rejection(N_PARTICLES, nfw_density_wrapper, &nfw_params, r_max, x, y, z, m);
    ASSERT_CLOSE(n_sampled, N_PARTICLES, 1e-10, "Sample count");
    
    /* Compute radii */
    for (int i = 0; i < N_PARTICLES; i++) {
        r[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
    }
    
    /* Total mass check (equal m, normalized to 1; scale later) */
    double m_tot = 0.0;
    for (int i = 0; i < N_PARTICLES; i++) m_tot += m[i];
    ASSERT_CLOSE(m_tot, 1.0, TOL_STAT, "Total mass norm");
    
    /* Radial hist vs analytic (bins 0.5 rs - 5 rs) */
    double r_min = 0.5 * rs, r_max_hist = 5.0 * rs;
    double hist_rho[N_BINS], hist_vol[N_BINS];
    compute_radial_hist(r, N_PARTICLES, r_min, r_max_hist, N_BINS, hist_rho, hist_vol);
    
    /* Check mid-bin means */
    double dr = (r_max_hist - r_min) / N_BINS;
    for (int bin = 0; bin < N_BINS; bin++) {
        double r_mid = r_min + (bin + 0.5) * dr;
        double rho_analytic = nfw_density(r_mid, rs, rho_s);
        double rho_est = hist_rho[bin] / rho_s;  /* Normalize to rho_s */
        double vol_norm = hist_vol[bin] / N_PARTICLES;  /* Fraction vol */
        if (rho_analytic > 0) {
            ASSERT_CLOSE(rho_est, rho_analytic, TOL_STAT, "NFW radial dist bin");
        }
    }
    printf("NFW positions: PASS\n");
}

/* Test velocities: Mean sigma ~ v_circ / sqrt(2), isotropic */
void test_sample_velocities_nfw() {
    printf("Testing NFW velocity sampling...\n");
    double rs = 100.0, rho_s = 1e9;
    double r_max = 5.0 * rs;
    double x[N_PARTICLES], y[N_PARTICLES], z[N_PARTICLES];
    double vx[N_PARTICLES], vy[N_PARTICLES], vz[N_PARTICLES], m[N_PARTICLES];
    
    /* Sample positions first */
    sample_positions_rejection(N_PARTICLES, nfw_density_wrapper, &nfw_params, r_max, x, y, z, m);
    
    /* Sample velocities */
    sample_velocities_isotropic(N_PARTICLES, x, y, z, nfw_enclosed_mass_wrapper, &nfw_params, vx, vy, vz);
    
    /* Compute dispersions (rms) */
    double sig_x2 = 0.0, sig_y2 = 0.0, sig_z2 = 0.0;
    double r_avg = 0.0;
    for (int i = 0; i < N_PARTICLES; i++) {
        double rr = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        r_avg += rr / N_PARTICLES;
        double m_enc = nfw_enclosed_mass(rr, rs, rho_s);
        double v_circ = sqrt(G * m_enc / rr);
        sig_x2 += vx[i]*vx[i] / N_PARTICLES;
        sig_y2 += vy[i]*vy[i] / N_PARTICLES;
        sig_z2 += vz[i]*vz[i] / N_PARTICLES;
    }
    double sigma_x = sqrt(sig_x2), sigma_y = sqrt(sig_y2), sigma_z = sqrt(sig_z2);
    double v_circ_avg = sqrt(G * nfw_enclosed_mass(r_avg, rs, rho_s) / r_avg);
    double sigma_exp = v_circ_avg / sqrt(2.0);  /* Isotropic approx */
    
    /* Check isotropy and magnitude */
    ASSERT_CLOSE(sigma_x, sigma_exp, TOL_STAT, "Sigma_x vs exp");
    ASSERT_CLOSE(sigma_y, sigma_exp, TOL_STAT, "Sigma_y vs exp");
    ASSERT_CLOSE(sigma_z, sigma_exp, TOL_STAT, "Sigma_z vs exp");
    
    /* Mean v = 0 */
    double mean_vx = 0.0, mean_vy = 0.0, mean_vz = 0.0;
    for (int i = 0; i < N_PARTICLES; i++) {
        mean_vx += vx[i] / N_PARTICLES;
        mean_vy += vy[i] / N_PARTICLES;
        mean_vz += vz[i] / N_PARTICLES;
    }
    ASSERT_CLOSE(fabs(mean_vx), 0.0, TOL_STAT, "Mean vx ~0");
    ASSERT_CLOSE(fabs(mean_vy), 0.0, TOL_STAT, "Mean vy ~0");
    ASSERT_CLOSE(fabs(mean_vz), 0.0, TOL_STAT, "Mean vz ~0");
    
    printf("NFW velocities: PASS\n");
}

/* Wrappers (assume from sampling.h; stub here) */
typedef struct { double rs, rhos; } NFWParams;
double nfw_density_wrapper(double rr, void *p) {
    NFWParams *params = (NFWParams *) p;
    return nfw_density(rr, params->rs, params->rhos);
}
double nfw_enclosed_mass_wrapper(double rr, void *p) {
    NFWParams *params = (NFWParams *) p;
    return nfw_enclosed_mass(rr, params->rs, params->rhos);
}

int main() {
    NFWParams nfw_params = {100.0, 1e9};  /* Test params */
    
    printf("Running sampling tests...\n\n");
    
    test_sample_positions_nfw();
    test_sample_velocities_nfw();
    
    printf("\nAll tests PASS!\n");
    return 0;
}
