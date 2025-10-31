# Architecture Overview

This document describes the high-level architecture of the Models-of-Galaxies-and-Clusters project, a modular Python framework for simulating and analyzing galaxy and cluster formation using semi-analytic models and N-body/hydrodynamic initial conditions. The project supports equilibrium models for dark matter halos, intracluster medium (ICM) gas, and stellar components, compatible with simulation codes like Gadget or Enzo. It emphasizes scalability for large-scale cosmological volumes and integration with observational data (e.g., from SDSS or DES). [web:2][web:46]

## Project Structure

The repository is organized into core directories for modularity, separation of concerns, and ease of extension:

- **models/**: Core simulation classes for physical components.
  - **galaxy.py**: Defines GalaxyModel for individual galaxies (disk, elliptical, spiral). Uses exponential or de Vaucouleurs profiles for surface brightness; integrates stellar dynamics via Jeans equations.
  - **cluster.py**: Implements ClusterModel combining Halo, Gas, and Stars. Supports NFW for dark matter and Beta models for gas density; normalizes to M200 virial mass.
  - **profiles.py**: Mathematical profiles (e.g., NFW, BetaModel) with methods for density ρ(r), enclosed mass M(<r), and potential Φ(r). Uses SciPy for numerical integration.
  - **halo.py**: Dark matter halo generator; particle sampling via rejection methods for NFW distribution.

- **scripts/**: Executable utilities for workflow automation.
  - **generate.py**: Main script to create initial conditions (ICs); generates particle positions/velocities and exports to HDF5 or binary formats.
  - **visualize.py**: Plotting tools using Matplotlib/Astropy; supports 2D projections, radial profiles, and density maps with logarithmic scaling.
  - **analyze.py**: Post-processing for scaling relations (e.g., L_X - T relation for clusters); computes baryon fractions and substructure via friends-of-friends algorithm.

- **data/**: Sample inputs and outputs.
  - Cosmology parameters (e.g., Planck18 h=0.674, Ω_m=0.315).
  - Pre-generated ICs for test clusters (e.g., M_vir=10^14-10^15 M_⊙ at z=0).

- **docs/**: Documentation and examples.
  - This ARCHITECTURE.md.
  - API references generated via Sphinx.

- **tests/**: Unit tests with pytest; covers profile normalizations and particle generation (e.g., assert enclosed_mass(r_200) ≈ M200).

- **config/**: YAML/JSON files for parameters (e.g., concentration c=10, β=0.67 for gas).

Root files: setup.py (pip installable), requirements.txt (NumPy, SciPy, Astropy, h5py, Matplotlib), and README.md with usage examples. [web:96][web:110]

## Data Flow

1. **Input**: User specifies cosmology (e.g., via config file), mass scales, and resolution (N_particles ~10^5-10^7).
2. **Model Initialization**: Cluster/GalaxyModel loads profiles; computes scale radii (r_s = R_vir / c) and normalizations (e.g., ρ_s from M200).
3. **Generation**: 
   - Halo particles sampled from ρ(r); velocities from Jeans or phase-space distribution.
   - Gas particles via Beta model; stars from King profiles for BCGs.
   - Combine into multi-type particle set (DM=1, gas=0, stars=4 in Gadget convention).
4. **Output**: HDF5 with positions, velocities, masses, types; or direct export to simulation ICs.
5. **Analysis/Visualization**: Scripts load outputs; compute observables (e.g., X-ray luminosity from gas pressure); plot projections or profiles.
6. **Extension Points**: Custom profiles via subclassing (e.g., add Einasto for halos); parallelization with MPI for large N.

Data types: Positions/velocities in kpc/km/s; masses in M_⊙; uses Astropy units for consistency. [web:123][web:124]

## Key Design Principles

- **Modularity**: Each component (halo, gas) is independent, allowing mix-and-match (e.g., NFW halo + isothermal gas).
- **Performance**: Vectorized NumPy operations; optional Cython for tight loops (e.g., potential evaluation). For large clusters, integrate with yt or SWIFTSIM for post-simulation viz.
- **Validation**: Profiles benchmarked against analytic solutions (e.g., M(<r) for NFW: 4π ρ_s r_s^3 [ln(1+x) - x/(1+x)] where x=r/r_s). Tested against observations like cluster baryon fractions f_bar ≈ 0.15. [web:46][web:50]
- **Scalability**: Supports hierarchical merging (semi-analytic) for cluster evolution; z-dependent via cosmology module (Astropy.cosmology).
- **Dependencies**:
  - Core: NumPy (>=1.21), SciPy (>=1.7), Astropy (>=5.0).
  - Viz: Matplotlib (>=3.5).
  - Optional: h5py for I/O, pytest for tests.

## Limitations and Future Work

Current version assumes spherical symmetry and hydrostatic equilibrium; no full hydrodynamics. Future: Add merger trees via semi-analytic galaxy formation (e.g., integrate with GALAXY_MAKER); support for weak lensing maps. Compare to observations like The300 project for 324 clusters. [web:2][web:46]

For contributions, see CONTRIBUTING.md. Contact: nayan993@example.com.