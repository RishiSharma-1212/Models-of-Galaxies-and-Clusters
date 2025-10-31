/**
 * test_integration.c
 *
 * Unit tests for integration.c: Trapezoidal, Simpson's, adaptive quadrature,
 * and enclosed mass numerical vs. analytic for NFW/Beta profiles.
 * Tests basic integrals (e.g., int x dx = x^2/2) and astronomy-specific
 * M(<r) = 4 pi int r^2 rho dr against known values.
 *
 * Tolerance: 1e-6 abs for basic, 1e-4 rel for mass (numerical error).
 * Uses assertions; prints PASS/FAIL, exits 1 on failure.
 *
 * References:
 * - Quadrature: Numerical Recipes (Press et al. 2007) [web:423][web:426]
 * - Enclosed mass tests: NFW analytic [web:424][web:425]
 * - Astronomy: Halo integration validation [web:427]
 *
 * Usage: gcc test_integration.c ../src/integration.o ../src/density.o -lm -o test_integration; ./test_integration
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/integration.h"
#include "../src/density.h"
#include "../src/utils.h"  /* For safe math */

/* Tolerance */
#define TOL_BASIC 1e-6
#define TOL_MASS 1e-4

/* Assertion macro */
#define ASSERT_CLOSE(val, exp, tol, msg) do { \
    double diff = fabs((val) - (exp)); \
    double rel = diff / (fabs(exp) + 1e-10); \
    if (rel > tol) { \
        fprintf(stderr, "FAIL %s: Got %.8g, exp %.8g (rel err=%.2e)\n", msg, val, exp, rel); \
        exit(1); \
    } \
} while(0)

/* Test integrand: f(x) = x (analytic int_a^b x dx = (b^2 - a^2)/2) */
double integrand_x(double x, void *params) {
    (void) params;  /* Unused */
    return x;
}

/* Test trapezoidal and Simpson's on simple integral */
void test_basic_quadrature() {
    printf("Testing basic quadrature...\n");
    double a = 0.0, b = 1.0;
    double exact = 0.5;  /* int_0^1 x dx */
    int N = 100;  /* Subintervals */
    
    /* Trapezoidal composite */
    double trap = trap_composite(a, b, N, integrand_x, NULL);
    ASSERT_CLOSE(trap, exact, TOL_BASIC, "Trapezoidal composite");
    
    /* Simpson composite (N even) */
    N = 100;  /* Even */
    double simp = simpson_composite(a, b, N, integrand_x, NULL);
    ASSERT_CLOSE(simp, exact, TOL_BASIC / 10.0, "Simpson composite");  /* Higher order */
    
    /* Single interval */
    double trap_single = trap_integral(a, b, integrand_x, NULL);
    ASSERT_CLOSE(trap_single, exact, TOL_BASIC, "Trapezoidal single");
    
    double simp_single = simpson_integral(a, b, integrand_x, NULL);
    ASSERT_CLOSE(simp_single, exact, TOL_BASIC, "Simpson single");
    
    printf("Basic quadrature: PASS\n");
}

/* Test adaptive Simpson */
void test_adaptive_simpson() {
    printf("Testing adaptive Simpson...\n");
    double a = 0.0, b = 2.0;
    double exact = 8.0 / 3.0;  /* int_0^2 x^2 dx */
    auto f_x2 = [](double x, void *p) { (void)p; return x * x; };
    double result;
    double abserr;
    int stat = adaptive_simpson(a, b, 1e-8, f_x2, NULL, &result, &abserr);
    ASSERT_CLOSE(stat, 0, 1e-10, "Adaptive status");  /* 0 success */
    ASSERT_CLOSE(result, exact, 1e-6, "Adaptive result");
    ASSERT_CLOSE(abserr, 0.0, 1e-6, "Adaptive error est < tol");
    printf("Adaptive Simpson: PASS\n");
}

/* Test numerical enclosed mass for NFW (vs analytic) */
void test_numerical_enclosed_mass_nfw() {
    printf("Testing numerical enclosed mass (NFW)...\n");
    double rs = 100.0, rho_s = 1e9;
    double r_test = rs;
    double m_analytic = nfw_enclosed_mass(r_test, rs, rho_s);
    double m_num = numerical_nfw_mass(r_test, rs, rho_s, 1000);  /* N=1000 */
    ASSERT_CLOSE(m_num, m_analytic, TOL_MASS, "NFW num vs analytic");
    
    /* At r=0 */
    m_num = numerical_nfw_mass(0.0, rs, rho_s, 1000);
    ASSERT_CLOSE(m_num, 0.0, TOL_MASS, "NFW num M(0)");
    
    /* Virial normalization */
    double rho_crit = compute_rho_crit(0.7);
    double m_vir = 1e15, c = 10.0;
    double rho_s_vir = compute_nfw_rho_s(m_vir, c, rho_crit);
    double r_vir = rs * c;  /* Assume rs = r_vir / c */
    m_analytic = nfw_enclosed_mass(r_vir, rs, rho_s_vir);
    m_num = numerical_nfw_mass(r_vir, rs, rho_s_vir, 1000);
    ASSERT_CLOSE(m_num, m_vir, 1e-3, "NFW virial num");  /* Allow 0.1% err */
    
    printf("NFW numerical mass: PASS\n");
}

/* Test Beta enclosed mass (no analytic, but consistency with low N increase) */
void test_numerical_enclosed_mass_beta() {
    printf("Testing numerical enclosed mass (Beta)...\n");
    double rc = 50.0, beta = 2.0/3.0, rho_0 = 1e10;
    double r_test = 2.0 * rc;
    double m_low = numerical_enclosed_mass(r_test, 100, beta_density_wrapper, &beta_params);  /* Stub params */
    double m_high = numerical_enclosed_mass(r_test, 1000, beta_density_wrapper, &beta_params);
    ASSERT_CLOSE(m_high, m_low, TOL_MASS / 10.0, "Beta convergence");  /* Higher N closer */
    
    printf("Beta numerical mass: PASS\n");
}

/* Wrapper for Beta (assume defined in integration.h or here) */
typedef struct { double rc, beta, rho0; } BetaParams;
double beta_density_wrapper(double r, void *p) {
    BetaParams *params = (BetaParams *) p;
    return beta_density(r, params->rc, params->beta, params->rho0);
}

int main() {
    BetaParams beta_params = {50.0, 2.0/3.0, 1e10};  /* Global for tests */
    
    printf("Running integration tests...\n\n");
    
    test_basic_quadrature();
    test_adaptive_simpson();
    test_numerical_enclosed_mass_nfw();
    test_numerical_enclosed_mass_beta();
    
    printf("\nAll tests PASS!\n");
    return 0;
}
