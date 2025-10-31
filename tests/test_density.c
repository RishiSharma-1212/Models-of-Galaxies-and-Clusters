/**
 * test_density.c
 *
 * Unit tests for density.c functions: NFW, Beta, Hernquist, Plummer profiles.
 * Verifies analytic forms at sample points, enclosed mass normalization,
 * and limits (r->0, r->inf). Uses tolerance 1e-6 for floating-point.
 *
 * Runs self-contained; prints PASS/FAIL and exits 0 on success.
 * Benchmarks: Known rho/M values from literature (e.g., NFW at r=rs).
 *
 * References:
 * - NFW: Navarro et al. (1996) [web:410][web:422]
 * - Hernquist: Hernquist (1990) [web:406]
 * - Plummer: Plummer (1911) [web:212]
 * - Beta: Cavaliere & Fusco-Femiano (1976) [web:157]
 *
 * Usage: gcc test_density.c ../src/density.o -lm -o test_density; ./test_density
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/density.h"  /* Adjust path as needed */
#include "../src/utils.h"   /* For safe math, rho_crit */

/* Tolerance for floating-point comparisons */
#define TOL 1e-6

/* Assertion macro: Print FAIL and exit if not close */
#define ASSERT_CLOSE(val, exp, tol, msg) do { \
    double diff = fabs((val) - (exp)); \
    if (diff > tol * fabs((exp) + 1e-10)) { \
        fprintf(stderr, "FAIL %s: Got %.8g, expected %.8g (diff=%.2e)\n", msg, val, exp, diff); \
        exit(1); \
    } \
} while(0)

/* Test NFW density */
void test_nfw_density() {
    printf("Testing NFW density...\n");
    double rs = 100.0;  /* kpc */
    double rho_s = 1e9; /* Msun/kpc^3 */
    
    /* At r=0: rho ~ rho_s / (0 * 1^2) -> inf, but impl clamps to 0? Test r small */
    ASSERT_CLOSE(nfw_density(1e-3 * rs, rs, rho_s), rho_s, TOL, "NFW r~0");
    
    /* At r=rs: rho = rho_s / (1 * 2^2) = rho_s / 4 */
    ASSERT_CLOSE(nfw_density(rs, rs, rho_s), rho_s / 4.0, TOL, "NFW at rs");
    
    /* At r=2 rs: x=2, rho = rho_s / (2 * 3^2) = rho_s / 18 */
    ASSERT_CLOSE(nfw_density(2.0 * rs, rs, rho_s), rho_s / 18.0, TOL, "NFW at 2rs");
    
    /* Large r: rho ~ rho_s * (rs/r)^3 -> small */
    ASSERT_CLOSE(nfw_density(1e3 * rs, rs, rho_s), 0.0, TOL, "NFW large r");
    printf("NFW density: PASS\n");
}

/* Test NFW enclosed mass */
void test_nfw_enclosed_mass() {
    printf("Testing NFW enclosed mass...\n");
    double rs = 100.0;
    double rho_s = 1e9;
    double m_rs = 4.0 * M_PI * rho_s * rs*rs*rs * (log(2.0) - 0.5);  /* Analytic at r=rs */
    
    ASSERT_CLOSE(nfw_enclosed_mass(rs, rs, rho_s), m_rs, TOL, "NFW M(<rs)");
    
    /* At r=0: M=0 */
    ASSERT_CLOSE(nfw_enclosed_mass(0.0, rs, rho_s), 0.0, TOL, "NFW M(0)");
    
    /* Virial check: Use compute_nfw_rho_s */
    double rho_crit = compute_rho_crit(0.7);
    double m_vir = 1e15;
    double c = 10.0;
    double rho_s_vir = compute_nfw_rho_s(m_vir, c, rho_crit);
    double r_vir = rs * c;  /* rs = r_vir / c */
    double m_vir_comp = nfw_enclosed_mass(r_vir, r_vir / c, rho_s_vir);
    ASSERT_CLOSE(m_vir_comp, m_vir, 1e-3, "NFW M_vir normalization");  /* Rel tol 0.1% */
    
    printf("NFW enclosed mass: PASS\n");
}

/* Test Hernquist density */
void test_hernquist_density() {
    printf("Testing Hernquist density...\n");
    double a = 50.0;  /* kpc */
    double M = 1e12;  /* Msun */
    
    /* At r=0: rho = M a / (2 pi r (r+a)^3) -> inf, but finite limit? Impl clamps */
    ASSERT_CLOSE(hernquist_density(1e-3 * a, a, M), (M * a) / (2.0 * M_PI * pow(a, 3.0) * 0.001), TOL, "Hernquist r~0 approx");
    
    /* At r=a: rho = M / (2 pi a^3) * (1/8) */
    ASSERT_CLOSE(hernquist_density(a, a, M), M / (16.0 * M_PI * a * a * a), TOL, "Hernquist at a");
    
    printf("Hernquist density: PASS\n");
}

/* Test Plummer density */
void test_plummer_density() {
    printf("Testing Plummer density...\n");
    double a = 100.0;  /* kpc */
    double M = 1e12;
    double rho0 = 3.0 * M / (4.0 * M_PI * pow(a, 3.0));  /* Central rho */
    
    /* At r=0: rho = 3 M / (4 pi a^3) */
    ASSERT_CLOSE(plummer_density(0.0, a, M), rho0, TOL, "Plummer at 0");
    
    /* At r=a: (1+1)^{-5/2} = 2^{-2.5} = 1/(4 sqrt(2)) */
    double factor = pow(2.0, -2.5);
    ASSERT_CLOSE(plummer_density(a, a, M), rho0 * factor, TOL, "Plummer at a");
    
    printf("Plummer density: PASS\n");
}

/* Test Beta density */
void test_beta_density() {
    printf("Testing Beta density...\n");
    double rc = 50.0;   /* kpc */
    double beta = 2.0/3.0;
    double rho0 = 1e10; /* Msun/kpc^3 */
    
    /* At r=0: rho = rho0 */
    ASSERT_CLOSE(beta_density(0.0, rc, beta, rho0), rho0, TOL, "Beta at 0");
    
    /* At r=rc: [1+1]^{-3 beta /2} = 2^{-1} = 0.5 for beta=2/3 */
    ASSERT_CLOSE(beta_density(rc, rc, beta, rho0), rho0 / 2.0, TOL, "Beta at rc");
    
    printf("Beta density: PASS\n");
}

/* Main test runner */
int main() {
    printf("Running density profile tests...\n\n");
    
    test_nfw_density();
    test_nfw_enclosed_mass();
    test_hernquist_density();
    test_plummer_density();
    test_beta_density();
    
    printf("\nAll tests PASS!\n");
    return 0;
}
