/**
 * integration.h
 *
 * Header file for integration.c: Numerical quadrature for astronomy profiles.
 * Declares trapezoidal, Simpson's, and adaptive methods for integrals like
 * enclosed mass M(<r) = 4 π ∫_0^r s^2 ρ(s) ds from density functions.
 *
 * Compatible with density.h for NFW/Beta/etc. Uses function pointers for generality.
 *
 * Units: r [kpc], integral [Msun] for mass.
 *
 * References:
 * - Numerical Recipes (Press et al. 2007) [web:298][web:305]
 * - QUADPACK/GSL for quadrature [web:301]
 * - Astronomy: Halo mass profiles [web:306][web:311]
 *
 * Usage: #include "integration.h"; #include "density.h"; compile with -lm.
 */

#ifndef INTEGRATION_H
#define INTEGRATION_H

/* Constants */
#define PI 3.14159265358979323846

/* Typedef for integrand function pointer (general: double f(double x, void *params)) */
typedef double (*IntegrandFunc)(double x, void *params);

/* Typedef for density function pointer (specific for mass integrands) */
typedef double (*DensityFunc)(double r, void *params);

/**
 * Trapezoidal rule prototypes
 */

/**
 * Single-interval trapezoidal: ∫_a^b f(x) dx ≈ (b-a)/2 * (f(a) + f(b))
 * @param a Lower limit
 * @param b Upper limit
 * @param f Integrand function
 * @param params Optional parameters
 * @return Integral approximation
 */
double trap_integral(double a, double b, IntegrandFunc f, void *params);

/**
 * Composite trapezoidal: With N subintervals (error O(h^2))
 * @param a Lower limit
 * @param b Upper limit
 * @param N Number of subintervals (>0)
 * @param f Integrand function
 * @param params Optional parameters
 * @return Integral approximation
 */
double trap_composite(double a, double b, int N, IntegrandFunc f, void *params);

/**
 * Simpson's rule prototypes
 */

/**
 * Single Simpson (parabola, N=2 subintervals): ∫_a^b f(x) dx ≈ (b-a)/6 * (f(a) + 4 f((a+b)/2) + f(b))
 * @param a Lower limit
 * @param b Upper limit
 * @param f Integrand function
 * @param params Optional parameters
 * @return Integral approximation
 */
double simpson_integral(double a, double b, IntegrandFunc f, void *params);

/**
 * Composite Simpson: With N even subintervals (error O(h^4))
 * @param a Lower limit
 * @param b Upper limit
 * @param N Number of subintervals (even, >0)
 * @param f Integrand function
 * @param params Optional parameters
 * @return Integral approximation
 */
double simpson_composite(double a, double b, int N, IntegrandFunc f, void *params);

/**
 * Mass-specific integrand: 4 π r^2 ρ(r), wrapper for density functions
 * @param r Radius [kpc]
 * @param params Profile parameters (e.g., rs, rho_s)
 * @param rho_func Density evaluator (e.g., nfw_density wrapper)
 * @return 4 π r^2 ρ(r) [Msun / kpc] for volume integral
 */
double mass_integrand(double r, void *params, DensityFunc rho_func);

/**
 * Numerical enclosed mass: M(<r_max) via composite Simpson on mass_integrand
 * @param r_max Upper radius [kpc]
 * @param N Subintervals (even, default 1000)
 * @param rho_func Density function (e.g., nfw_density with params)
 * @param params Profile struct (e.g., {rs, rho_s})
 * @return M(<r_max) [Msun]
 */
double numerical_enclosed_mass(double r_max, int N, DensityFunc rho_func, void *params);

/**
 * Example: Numerical NFW mass (wrapper for validation)
 * @param r Radius [kpc]
 * @param rs Scale [kpc]
 * @param rho_s Char density [Msun/kpc^3]
 * @param N Subintervals
 * @return M(<r) [Msun]
 */
double numerical_nfw_mass(double r, double rs, double rho_s, int N);

/* Adaptive Simpson prototypes (recursive for error control) */

/**
 * Adaptive Simpson quadrature: Subdivides until error < tol
 * @param a Lower limit
 * @param b Upper limit
 * @param tol Tolerance (absolute error)
 * @param f Integrand function
 * @param params Optional parameters
 * @param *result Output: Computed integral
 * @param *abserr Output: Estimated error
 * @return 0 on success, 1 on failure (e.g., max depth)
 */
int adaptive_simpson(double a, double b, double tol, IntegrandFunc f, void *params, double *result, double *abserr);

/* Additional helpers (e.g., for Beta or custom profiles) */
/* Similar prototypes for numerical_beta_mass, etc. */

#endif /* INTEGRATION_H */
