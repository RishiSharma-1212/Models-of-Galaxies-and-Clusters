/**
 * sampling.c
 *
 * Particle sampling for density profiles in galaxy/cluster N-body initial conditions.
 * Uses rejection sampling to generate positions according to rho(r) (e.g., NFW, Hernquist).
 * Draws uniform r, theta, phi; accepts if rho(r_sample) / rho_max > u ~ Uniform[0,1].
 * Velocities via Jeans equation stub (isotropic, sigma from Phi); extend for full DF.
 *
 * Integrates density.h; assumes spherical symmetry.
 * Output: Cartesian positions [kpc], masses [Msun]; redirect to file for ICs.
 *
 * Sample: NFW halo with N=1e5 particles, M_tot normalized.
 *
 * References:
 * - Rejection sampling: Press et al. (1992) [web:356]
 * - IC generation: Sirko (2005); Joyce et al. (2005) [web:357][web:366]
 * - N-body ICs: L'Huillier et al. (2014) [web:315]
 *
 * Usage: #include "sampling.h"; link with density.o -lm.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "density.h"

/* Constants */
#define PI 3.14159265358979323846
#define G 4.302e-3  /* km^2 s^{-2} kpc Msun^{-1} */

/* RNG seed */
static unsigned int seed = 0;

/**
 * Simple RNG: Park-Miller (minimal standard); replace with GSL for production.
 * @return Uniform [0,1)
 */
double rand_uniform() {
    const long int a = 16807, m = 2147483647, q = 127773, r = 2836;
    long int hi = seed / q;
    long int lo = seed % q;
    long int test = a * lo - r * hi;
    if (test > 0) {
        seed = test;
    } else {
        seed = test + m;
    }
    return (double) seed / m;
}

/**
 * Rejection sampling for positions in spherical shell; uniform in volume, accept w/ prob rho(r)/rho_max.
 * @param N_particles Number to sample
 * @param rho_func Density function (e.g., nfw_density)
 * @param params Profile params (rs, rho_s)
 * @param r_max Max radius [kpc] (enclosing volume)
 * @param *x, *y, *z Output position arrays [kpc]
 * @param *m Output mass array [Msun] (equal mass)
 * @return Number sampled (==N if success)
 */
int sample_positions_rejection(int N, double (*rho_func)(double, void*), void *params,
                               double r_max, double *x, double *y, double *z, double *m) {
    if (N <= 0 || r_max <= 0.0) return 0;

    /* Compute rho_max (brute force or analytic; here grid search for generality) */
    int n_grid = 1000;
    double dr = r_max / n_grid;
    double rho_max = 0.0;
    for (int i = 1; i < n_grid; i++) {
        double rr = i * dr;
        double rho = rho_func(rr, params);
        if (rho > rho_max) rho_max = rho;
    }
    if (rho_max <= 0.0) return 0;

    /* Enclosing volume V = 4/3 pi r_max^3; prob density 3/(r_max^3) for r^3 uniform */
    double vol = 4.0 / 3.0 * PI * r_max * r_max * r_max;
    double m_part = 1.0;  /* Normalize total M=1; scale later */

    int n_sampled = 0;
    while (n_sampled < N) {
        /* Propose uniform in volume: r^3 ~ U[0, V], theta~U[0,pi], phi~U[0,2pi] */
        double u_r3 = rand_uniform();  /* U[0,1) */
        double r_prop = r_max * cbrt(u_r3);
        double theta = acos(2.0 * rand_uniform() - 1.0);  /* Cos theta uniform */
        double phi = 2.0 * PI * rand_uniform();

        /* Evaluate rho at proposed r */
        double rho_prop = rho_func(r_prop, params);

        /* Acceptance prob = rho_prop / (rho_max * (3 / r_max^3)) but since uniform vol, simplify to rho_prop / rho_max * normalization */
        double u_accept = rand_uniform();
        if (u_accept < rho_prop / rho_max) {
            /* Accept: Convert to Cartesian */
            double sin_theta = sin(theta);
            x[n_sampled] = r_prop * sin_theta * cos(phi);
            y[n_sampled] = r_prop * sin_theta * sin(phi);
            z[n_sampled] = r_prop * cos(theta);
            m[n_sampled] = m_part;  /* Equal mass particles */

            n_sampled++;
        }
    }

    return n_sampled;
}

/**
 * Sample velocities (isotropic Jeans approximation): sigma^2(r) = 1/rho int_r^inf rho dPhi/dr dr (simplified).
 * Stub: Use circular velocity v_circ = sqrt(G M(<r)/r); Gaussian sigma = v_circ / sqrt(2) isotropic.
 * @param N_particles Number
 * @param x,y,z Positions [kpc]
 * @param rho_func Density func
 * @param mass_func Enclosed mass func (e.g., nfw_enclosed_mass)
 * @param params Profile params
 * @param *vx, *vy, *vz Output velocity arrays [km/s]
 * @return Void
 */
void sample_velocities_isotropic(int N, double *x, double *y, double *z,
                                 double (*mass_func)(double, void*),
                                 void *params, double *vx, double *vy, double *vz) {
    for (int i = 0; i < N; i++) {
        double rx = x[i], ry = y[i], rz = z[i];
        double r = sqrt(rx*rx + ry*ry + rz*rz);
        if (r <= 0.0) continue;

        double m_enc = mass_func(r, params);
        double v_circ = sqrt(G * m_enc / r);  /* km/s */
        double sigma = v_circ / sqrt(2.0);  /* Approx isotropic dispersion */

        /* Gaussian velocities: mean 0, sigma in each dir (isotropic) */
        double u_x = rand_uniform(), u_y = rand_uniform(), u_z = rand_uniform();
        vx[i] = sigma * sqrt(-2.0 * log(1.0 - u_x));  /* Box-Muller stub; full impl needed */
        vy[i] = sigma * sqrt(-2.0 * log(1.0 - u_y));
        vz[i] = sigma * sqrt(-2.0 * log(1.0 - u_z));
        /* Randomize signs/directions */
        if (rand_uniform() > 0.5) vx[i] = -vx[i];
        if (rand_uniform() > 0.5) vy[i] = -vy[i];
        if (rand_uniform() > 0.5) vz[i] = -vz[i];
    }
}

/**
 * Wrapper for NFW particle sampling (positions + velocities).
 * @param N_particles
 * @param rs Scale [kpc]
 * @param rho_s Char density [Msun/kpc^3]
 * @param r_max
 * @param *pos, *vel Output (3N arrays or struct)
 * @param *m Masses [Msun]
 * @return Success (1)
 */
int sample_nfw_particles(int N, double rs, double rho_s, double r_max,
                         double *x, double *y, double *z, double *vx, double *vy, double *vz, double *m) {
    typedef struct { double rs, rhos; } NFWParams;
    NFWParams p = {rs, rho_s};
    auto nfw_rho = [](double r, void *pp) {
        NFWParams *p = (NFWParams*) pp;
        return nfw_density(r, p->rs, p->rhos);
    };
    auto nfw_mass = [](double r, void *pp) {
        NFWParams *p = (NFWParams*) pp;
        return nfw_enclosed_mass(r, p->rs, p->rhos);
    };

    srand(time(NULL));  /* Seed RNG */
    seed = (unsigned int) time(NULL);

    int n_pos = sample_positions_rejection(N, nfw_rho, &p, r_max, x, y, z, m);
    if (n_pos == N) {
        sample_velocities_isotropic(N, x, y, z, nfw_mass, &p, vx, vy, vz);
        return 1;
    }
    return 0;
}

/* Similar wrappers for Hernquist, Beta gas... e.g., sample_hernquist_particles */

/**
 * Output particles to CSV/binary (stub: print positions/masses; extend for Gadget format).
 * @param N
 * @param x,y,z,vx,vy,vz,m Arrays
 */
void output_particles_csv(int N, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *m) {
    printf("x_kpc,y_kpc,z_kpc,vx_kms,vy_kms,vz_kms,m_Msun\n");
    for (int i = 0; i < N; i++) {
        printf("%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
               x[i], y[i], z[i], vx[i], vy[i], vz[i], m[i]);
    }
}

/* Main demo stub (call from main.c) */
/* int main() { ... sample_nfw_particles(10000, rs, rho_s, r_max, ...); output... } */
