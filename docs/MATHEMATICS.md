# Mathematical Foundations

This document details the core mathematical models used in the Models-of-Galaxies-and-Clusters project. We focus on spherical symmetry and hydrostatic equilibrium assumptions, common in semi-analytic and initial condition generation for N-body/hydrodynamic simulations. All models are implemented with Astropy units for dimensional consistency. Derivations are grounded in ΛCDM cosmology, where dark matter halos form hierarchically from initial Gaussian fluctuations. [web:158][web:165]

## Dark Matter Halos: Navarro-Frenk-White (NFW) Profile

The NFW profile describes the density distribution in dark matter halos, derived from high-resolution N-body simulations. It assumes a cuspy inner slope and shallower outer fall-off, universal across masses and redshifts. [web:153][web:159]

### Density Profile

The three-dimensional density ρ(r) is given by:

\[ \rho(r) = \frac{\rho_s}{(r / r_s) (1 + r / r_s)^2} \]

where:
- \( r_s \): Scale radius (kpc), related to virial radius \( R_{vir} \) by \( r_s = R_{vir} / c \) with concentration \( c \approx 5-20 \) (mass-dependent).
- \( \rho_s \): Characteristic density (M_⊙ kpc^{-3}), set by normalization.

The inner cusp (\( r \ll r_s \)) behaves as \( \rho \propto r^{-1} \), while outer (\( r \gg r_s \)) as \( \rho \propto r^{-3} \). This fits observations of rotation curves and lensing in clusters. [web:153][web:165]

### Normalization and Enclosed Mass

The characteristic density \( \rho_s \) is normalized such that the mass within the virial radius \( M_{vir} \) matches the halo mass:

\[ \rho_s = \frac{200 \rho_{crit}}{3 c^3 \left[ \ln(1 + c) - \frac{c}{1 + c} \right]} \]

with critical density \( \rho_{crit} = \frac{3 H(z)^2}{8 \pi G} \), H(z) from cosmology (e.g., Planck18). The enclosed mass M(<r) is:

\[ M(<r) = 4 \pi \rho_s r_s^3 \left[ \ln\left(1 + \frac{r}{r_s}\right) - \frac{r / r_s}{1 + r / r_s} \right] \]

At \( r = R_{vir} \), M(<R_{vir}) = M_{vir} \). For numerical stability, integrate via SciPy.quad. [web:159][web:168]

### Gravitational Potential

The potential Φ(r) (up to -G M / r) for NFW is:

\[ \Phi(r) = -4 \pi G \rho_s r_s^2 \frac{\ln(1 + r / r_s)}{r / r_s} \]

Used for velocity dispersion in Jeans equation for equilibrium. [web:153]

## Intracluster Medium: Beta Model

The Beta model parameterizes the ICM gas density, assuming isothermal spheres in hydrostatic equilibrium with the NFW potential. It fits X-ray surface brightness observations well. [web:154][web:156][web:163]

### Density Profile

\[ \rho_{gas}(r) = \rho_0 \left[ 1 + \left( \frac{r}{r_c} \right)^2 \right]^{-3 \beta / 2} \]

where:
- \( r_c \): Core radius (kpc, ~10-100 kpc).
- \( \beta \approx 2/3 \): Slope parameter, linking gas temperature to halo potential (\( \beta = \mu m_p \sigma^2 / k T \), σ velocity dispersion).
- \( \rho_0 \): Central density (M_⊙ kpc^{-3}), normalized to gas mass fraction f_gas ≈ 0.15.

The model has a flat core (r << r_c) and \( \rho \propto r^{-3\beta} \) asymptotically (~r^{-2} for β=2/3). [web:154][web:156]

### Normalization

Total gas mass M_gas = f_gas M_vir, with enclosed gas mass:

\[ M_{gas}(<r) = 4 \pi \rho_0 r_c^3 \left(1 - \frac{1}{(3\beta - 1/2)} \, _2F_1\left(\frac{3\beta}{2}, \frac{1}{2}; \frac{3\beta + 1/2}{2}; -\left(\frac{r}{r_c}\right)^2 \right) \right) \]

For numerical use, integrate \( 4\pi r^2 \rho(r) \). Hypergeometric \( _2F_1 \) for analytic form. [web:156]

### Pressure and Temperature

Assuming isothermal gas at T (~10^7-10^8 K), pressure P(r):

\[ P(r) = \frac{\rho_{gas}(r) k T}{\mu m_p} \]

with μ ≈ 0.59, k Boltzmann constant, m_p proton mass. Hydrostatic equilibrium: dP/dr = -ρ_gas g(r), g(r) = -dΦ/dr from NFW. [web:162]

## Galaxy Models: Hernquist and Exponential Disks

For central galaxies (BCGs), use Hernquist for spheroids (deprojected de Vaucouleurs):

\[ \rho(r) = \frac{M a}{2\pi r (r + a)^3} \]

a scale parameter. For disks: surface density Σ(R) = Σ_0 exp(-R / R_d), deprojected via integrals. [web:157]

## Hydrostatic Equilibrium and Virial Theorem

Cluster equilibrium: ∇P = -ρ ∇Φ, integrated for M(<r):

\[ M(<r) = -\frac{r^2}{\rho_{gas} G} \frac{d \rho_{gas}}{dr} \frac{k T}{\mu m_p} \]

Combines X-ray data (ρ_gas from emissivity ∝ ρ^2) with lensing (total M). [web:162][web:170]

## Cosmological Context

Halo mass function from Press-Schechter or excursion sets; concentration c(M,z) ≈ 7.85 (M / 10^{14} M_⊙)^{-0.081} (1+z)^{0.71}. Growth via D(z) in simulations. [web:158][web:165]

For derivations, see Navarro et al. (1996) [web:165] and Cavaliere & Fusco-Femiano (1976) [web:154]. Implementations in profiles.py use these for vectorized computation.