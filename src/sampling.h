/**
 * sampling.h
 *
 * Header file for sampling.c: Particle initial conditions via rejection sampling
 * for density profiles in N-body galaxy/cluster simulations.
 * Declares functions for positions (rejection method) and velocities (isotropic Jeans approx),
 * supporting NFW, Hernquist, etc., from density.h.
 *
 * Assumes spherical symmetry; outputs positions/velocities for ICs.
 * RNG: Simple uniform; use GSL for production.
 *
 * Units: Positions [kpc], velocities [km/s], masses [Msun].
 *
 * References:
 * - Rejection sampling: Von Neumann (1951); Press et al. (1992) [web:356]
 * - N-body ICs: Sirko (2005); Angulo et al. (2013) [web:357][web:361]
 * - Applications: Halo sampling in simulations [web:353][web:366]
 *
 * Usage: #include "sampling.h"; #include "density.h"; link -lm.
 */

#ifndef SAMPLING_H
#define SAMPLING_H

#include "density.h"  /* For density/mass functions */

/* Constants */
#define PI 3.14159265358979323846
#define G 4.302e-3  /* km^2 s^{-2} kpc Msun^{-1} */

/* Typedef for density and enclosed mass function pointers */
typedef double (*DensityFunc)(double r, void *params);
typedef double (*EnclosedMassFunc)(double r, void *params);

/**
 * Simple RNG prototypes (Park-Miller; seed externally)
 */

/**
 * Uniform random [0,1)
 * @return u ~ U[0,1)
 */
double rand_uniform(void);

/**
 * Rejection sampling for positions: Propose uniform in sphere of r_max, accept w/ p = rho(r)/rho_max
 * @param N_particles Number to sample
 * @param rho_func Density evaluator
 * @param params Profile parameters (e.g., rs, rho_s)
 * @param r_max Enclosing radius [kpc]
 * @param *x, *y, *z Output positions [kpc] (pre-allocated arrays of size N)
 * @param *m Output masses [Msun] (equal mass, total M=1 normalized)
 * @return Number successfully sampled (==N if no failure)
 */
int sample_positions_rejection(int N, DensityFunc rho_func, void *params,
                               double r_max, double *x, double *y, double *z, double *m);

/**
 * Isotropic velocity sampling (approximate Jeans): Gaussian sigma ~ v_circ / sqrt(2)
 * @param N_particles Number
 * @param x,y,z Positions [kpc]
 * @param mass_func Enclosed mass evaluator
 * @param params Profile parameters
 * @param *vx, *vy, *vz Output velocities [km/s] (pre-allocated arrays)
 * @return Void
 */
void sample_velocities_isotropic(int N, double *x, double *y, double *z,
                                 EnclosedMassFunc mass_func, void *params,
                                 double *vx, double *vy, double *vz);

/**
 * NFW-specific particle sampler: Positions + velocities
 * @param N_particles
 * @param rs Scale radius [kpc]
 * @param rho_s Characteristic density [Msun/kpc^3]
 * @param r_max Enclosing radius [kpc]
 * @param *x,y,z Positions [kpc]
 * @param *vx,vy,vz Velocities [km/s]
 * @param *m Masses [Msun]
 * @return 1 on success, 0 on failure
 */
int sample_nfw_particles(int N, double rs, double rho_s, double r_max,
                         double *x, double *y, double *z,
                         double *vx, double *vy, double *vz, double *m);

/**
 * Output particles to CSV (positions, velocities, masses)
 * @param N Number
 * @param x,y,z,vx,vy,vz,m Arrays
 * @return Void (prints to stdout)
 */
void output_particles_csv(int N, double *x, double *y, double *z,
                          double *vx, double *vy, double *vz, double *m);

/* Optional: Wrappers for other profiles, e.g., sample_hernquist_particles */
/* int sample_hernquist_particles(int N, double a, double M, double r_max, ...); */

/* RNG seeding */
/**
 * Set RNG seed
 * @param s Seed value (e.g., time(NULL))
 * @return Void
 */
void srand_custom(unsigned int s);

#endif /* SAMPLING_H */
