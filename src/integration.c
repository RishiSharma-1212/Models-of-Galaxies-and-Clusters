/**
 * integration.c
 *
 * Numerical integration utilities for computing enclosed mass and other quantities
 * from density profiles in galaxy/cluster models. Uses trapezoidal and Simpson's rules
 * for radial integrals, e.g., M(<r) = 4 π ∫_0^r s^2 ρ(s) ds.
 *
 * Integrates with density.h for profile evaluations. Suitable for NFW, Beta, etc.,
 * where analytic M(<r) may be approximate or unavailable.
 *
 * Units: r [kpc], rho [Msun/kpc^3], M [Msun].
 *
 * References:
 * - Numerical Recipes (Press et al. 2007) for quadrature methods [web:282]
 * - Astronomy applications: Enclosed mass in halos [web:279][web:280]
 *
 * Usage: #include "integration.h"; link with density.o -lm.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "density.h"
#include "integration.h"  // Assumed header with prototypes

/* Constants */
#define PI 3.14159265358979323846

/**
 * Trapezoidal rule for ∫_a^b f(x) dx, single interval.
 * @param a Lower limit
 * @param b Upper limit
 * @param f Function pointer (e.g., integrand)
 * @return Integral value
 */
double trap_integral(double a, double b, double (*f)(double)) {
    if (a >= b) return 0.0;
    return (b - a) * 0.5 * (f(a) + f(b));
}

/**
 * Composite trapezoidal rule for ∫_a^b f(x) dx with N subintervals.
 * @param a Lower limit
 * @param b Upper limit
 * @param N Number of subintervals (even for accuracy)
 * @param f Function pointer
 * @return Integral value
 */
double trap_composite(double a, double b, int N, double (*f)(double)) {
    if (N <= 0 || a >= b) return 0.0;
    double h = (b - a) / N;
    double sum = 0.5 * (f(a) + f(b));
    for (int i = 1; i < N; i++) {
        sum += f(a + i * h);
    }
    return h * sum;
}

/**
 * Simpson's rule for ∫_a^b f(x) dx, single parabola (N=2 subintervals).
 * @param a Lower limit
 * @param b Upper limit
 * @param f Function pointer
 * @return Integral value
 */
double simpson_integral(double a, double b, double (*f)(double)) {
    if (a >= b) return 0.0;
    double c = (a + b) / 2.0;
    return (b - a) / 6.0 * (f(a) + 4.0 * f(c) + f(b));
}

/**
 * Composite Simpson's rule for ∫_a^b f(x) dx with N subintervals (N even).
 * @param a Lower limit
 * @param b Upper limit
 * @param N Number of subintervals (must be even)
 * @param f Function pointer
 * @return Integral value
 */
double simpson_composite(double a, double b, int N, double (*f)(double)) {
    if (N % 2 != 0 || N <= 0 || a >= b) return 0.0;
    double h = (b - a) / N;
    double sum = f(a) + f(b);
    double sum_odd = 0.0, sum_even = 0.0;
    for (int i = 1; i < N; i++) {
        double x = a + i * h;
        if (i % 2 == 0) {
            sum_even += f(x);
        } else {
            sum_odd += f(x);
        }
    }
    return h / 3.0 * (sum + 4.0 * sum_odd + 2.0 * sum_even);
}

/**
 * Integrand for enclosed mass: 4 π r^2 ρ(r), where ρ(r) from profile func.
 * Wrapper: Calls density func (e.g., nfw_density) to get ρ(r).
 * @param r Radius [kpc]
 * @param rho_func Pointer to density function (e.g., double (*)(double, double, double))
 * @param params Profile params (rs, rho_s, etc.; cast to void* for generality)
 * @return 4 π r^2 ρ(r) [Msun / kpc]
 */
double mass_integrand(double r, void *params, double (*rho_func)(double, void*)) {
    if (r <= 0.0) return 0.0;
    double rho = rho_func(r, params);
    return 4.0 * PI * r * r * rho;
}

/**
 * Compute enclosed mass M(<r_max) via numerical integration of mass_integrand.
 * Uses composite Simpson (preferred for smooth profiles).
 * @param r_max Upper radius [kpc]
 * @param N Subintervals (even, default 1000 for accuracy ~1e-4)
 * @param rho_func Density function pointer
 * @param params Profile params
 * @return M(<r_max) [Msun]
 */
double numerical_enclosed_mass(double r_max, int N, double (*rho_func)(double, void*), void *params) {
    if (r_max <= 0.0 || N <= 0) return 0.0;
    // Integrand wrapper: lambda-like via pointer
    auto integrand = [&](double r) { return mass_integrand(r, params, rho_func); };
    // For C, define static func or use callback; here assume rho_func takes void*
    return simpson_composite(0.0, r_max, N, integrand);  // Note: In pure C, replace lambda with explicit func
    // Fallback: Trapezoidal if N odd
    // return trap_composite(0.0, r_max, N, integrand);
}

/**
 * Example wrapper for NFW enclosed mass numerical (for validation vs analytic).
 * @param r Radius [kpc]
 * @param rs Scale [kpc]
 * @param rho_s Char density [Msun/kpc^3]
 * @param N Subintervals
 * @return M(<r) [Msun]
 */
double numerical_nfw_mass(double r, double rs, double rho_s, int N) {
    typedef struct { double rs, rho_s; } NFWParams;
    NFWParams p = {rs, rho_s};
    auto nfw_rho = [](double rr, void *pptr) {
        NFWParams *p = (NFWParams*) pptr;
        return nfw_density(rr, p->rs, p->rho_s);
    };
    return numerical_enclosed_mass(r, N, nfw_rho, &p);
}

/* Similar wrappers for Beta, Hernquist, Plummer... (e.g., numerical_beta_mass) */

/**
 * Adaptive quadrature stub (simple bisection for error control; for production, use full adaptive like from NR).
 * Recursively subdivides if error > tol.
 * @param a Lower
 * @param b Upper
 * @param tol Tolerance
 * @param f Func
 * @param *result Output integral
 * @param *abserr Output error estimate
 * @return Status (0 success)
 */
int adaptive_simpson(double a, double b, double tol, double (*f)(double), double *result, double *abserr) {
    double c = (a + b) / 2.0;
    double left = simpson_integral(a, c, f);
    double right = simpson_integral(c, b, f);
    double full = simpson_integral(a, b, f);
    *abserr = fabs((left + right) - full) / 15.0;  // Simpson error est.
    if (*abserr < tol) {
        *result = left + right;
        return 0;
    }
    double err_left, err_right;
    double res_left, res_right;
    int stat_left = adaptive_simpson(a, c, tol / 2.0, f, &res_left, &err_left);
    int stat_right = adaptive_simpson(c, b, tol / 2.0, f, &res_right, &err_right);
    if (stat_left == 0 && stat_right == 0) {
        *result = res_left + res_right;
        *abserr = err_left + err_right;
        return 0;
    }
    return 1;  // Failure
}

/* Vectorized integration for array of r_max (e.g., cumulative M(r)) */
/* Implement loop over numerical_enclosed_mass for simplicity */

