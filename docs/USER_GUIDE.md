# User Guide

This guide walks through installing, running, and extending the Models-of-Galaxies-and-Clusters project, a Python-based toolkit for generating equilibrium models of galaxies and clusters suitable for N-body and hydrodynamic simulations. It supports creating initial conditions (ICs) with dark matter halos, intracluster gas, and stellar components, exportable to formats like HDF5 or Gadget binary. [web:6][web:173] The project assumes basic familiarity with Python and astronomy concepts like NFW profiles; see MATHEMATICS.md for theory. [web:159][web:165]

## Installation

### Prerequisites

- Python 3.8+.
- NumPy, SciPy, Astropy (for units/cosmology), h5py (for HDF5 I/O), Matplotlib (for visualization).
- Optional: Gadget-2/3 or Enzo for running exported ICs. [web:6][web:110]

### From Source

1. Clone the repository:
