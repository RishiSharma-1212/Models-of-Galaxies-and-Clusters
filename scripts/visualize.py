"""
Visualization script for Models-of-Galaxies-and-Clusters project.
Plots 2D projections of galaxy and cluster models, including density profiles,
stellar distributions, and gas maps. Supports saving figures to PNG/PDF.

Usage:
    python visualize.py --model galaxy --input data/galaxy_sim.h5 --output plots/
    python visualize.py --model cluster --radius 500 --output plots/

Dependencies: numpy, matplotlib, astropy, h5py (for loading simulations).
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy import units as u
from astropy.coordinates import SkyCoord
import h5py  # For loading HDF5 simulation outputs
from models.galaxy import GalaxyModel  # Assume imported from project models/galaxy.py
from models.cluster import ClusterModel   # Assume imported from project models/cluster.py
from models.profiles import NFW, BetaModel  # For profile plotting

# Set global style for plots
plt.style.use('dark_background')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12


def load_simulation_data(filename):
    """
    Load particle positions, velocities, and masses from HDF5 simulation output.
    :param filename: Path to HDF5 file (e.g., 'data/galaxy_sim.h5')
    :return: dict with 'positions', 'masses', 'types' (0=gas, 1=DM, 2=stars)
    """
    data = {}
    with h5py.File(filename, 'r') as f:
        data['positions'] = np.array(f['positions'])  # Shape: (N_particles, 3)
        data['masses'] = np.array(f['masses'])        # Shape: (N_particles,)
        data['types'] = np.array(f['types'])          # Shape: (N_particles,)
        data['velocities'] = np.array(f['velocities']) if 'velocities' in f else None
    return data


def project_to_2d(positions, projection_axis=2, extent=(-500, 500, -500, 500)):
    """
    Project 3D positions to 2D sky plane (e.g., RA-Dec or x-y).
    :param positions: (N, 3) array in kpc
    :param projection_axis: Axis to project along (default z=2)
    :param extent: Plot extent in kpc
    :return: 2D positions and density map
    """
    if projection_axis == 2:
        proj_pos = positions[:, :2]
    else:
        # General projection (simplified)
        proj_pos = np.column_stack((positions[:, 0], positions[:, 1]))
    
    # Bin to density map
    hist, xedges, yedges = np.histogram2d(proj_pos[:, 0], proj_pos[:, 1],
                                          bins=256, range=[extent[:2], extent[2:]],
                                          weights=None)  # Density without weights for count
    return proj_pos, hist.T  # Transpose for imshow


def plot_galaxy(model_type='disk', radius=50, output_dir='plots/', save=True):
    """
    Visualize a single galaxy model (e.g., disk or elliptical).
    :param model_type: 'disk', 'spiral', or 'elliptical'
    :param radius: Galaxy radius in kpc
    :param output_dir: Directory to save plots
    :param save: Whether to save the figure
    """
    # Initialize model
    galaxy = GalaxyModel(type=model_type, radius=radius)
    
    # Generate 2D surface density profile (exponential disk example)
    r = np.linspace(0, radius, 1000)
    if model_type == 'disk':
        sigma = galaxy.surface_density(r)  # Assume method returns Sigma(r)
    else:
        sigma = galaxy.profile.density(r * u.kpc) * (r * u.kpc).to(u.cm, equivalencies=u.dimensional_equivalencies())  # Placeholder
    
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(r, sigma, 'w-', linewidth=2, label=f'{model_type.capitalize()} Galaxy')
    ax.set_rmax(radius)
    ax.set_rticks([10, 20, 30, 40, 50])
    ax.set_title(f'{model_type.capitalize()} Galaxy Surface Density Profile', color='white')
    ax.grid(True)
    
    if save:
        plt.savefig(f'{output_dir}/galaxy_{model_type}.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'{output_dir}/galaxy_{model_type}.pdf', bbox_inches='tight')
    plt.show()


def plot_cluster_radius(radius=500, profile='NFW', output_dir='plots/', save=True):
    """
    Plot radial profiles for a cluster (density, temperature, etc.).
    :param radius: Cluster virial radius in kpc
    :param profile: 'NFW' for DM or 'beta' for gas
    :param output_dir: Directory to save plots
    :param save: Whether to save the figure
    """
    r = np.linspace(0.1, radius, 1000) * u.kpc
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    if profile == 'NFW':
        prof = NFW(rs=100 * u.kpc, c=10)
        density = prof.density(r)
        enclosed_mass = prof.enclosed_mass(r)
        ax1.loglog(r.value, density.value, 'r-', label='DM Density')
        ax1.set_xlabel('Radius (kpc)')
        ax1.set_ylabel(r'Density (M$_\odot$ kpc$^{-3}$)')
        ax1.legend()
        
        ax2.semilogx(r.value, enclosed_mass.value, 'g-', label='Enclosed Mass')
        ax2.set_xlabel('Radius (kpc)')
        ax2.set_ylabel(r'Mass (M$_\odot$)')
        ax2.legend()
    elif profile == 'beta':
        gas_prof = BetaModel(rc=50 * u.kpc, beta=0.67)
        rho_gas = gas_prof.density(r)
        # Placeholder temperature profile (isothermal)
        T = np.full_like(r, 8e7)  # K
        pressure = gas_prof.pressure(r, T=T)
        
        ax1.loglog(r.value, rho_gas.value, 'b-', label='Gas Density')
        ax1.set_xlabel('Radius (kpc)')
        ax1.set_ylabel(r'Density (M$_\odot$ kpc$^{-3}$)')
        ax1.legend()
        
        ax2.loglog(r.value, pressure.value, 'm-', label='Gas Pressure')
        ax2.set_xlabel('Radius (kpc)')
        ax2.set_ylabel(r'Pressure (dyn cm$^{-2}$)')
        ax2.legend()
    
    plt.suptitle(f'Cluster {profile.upper()} Profile (R_vir = {radius} kpc)', color='white')
    
    if save:
        plt.savefig(f'{output_dir}/cluster_{profile}_profiles.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'{output_dir}/cluster_{profile}_profiles.pdf', bbox_inches='tight')
    plt.show()


def plot_simulation_projection(data_file, particle_type=1, extent=(-200, 200, -200, 200),
                               output_dir='plots/', save=True):
    """
    Visualize 2D projection of simulation particles (e.g., DM or gas map).
    :param data_file: HDF5 file with simulation data
    :param particle_type: Filter by type (0=gas, 1=DM, 2=stars)
    :param extent: Plot extent (kpc)
    :param output_dir: Directory to save plots
    :param save: Whether to save the figure
    """
    data = load_simulation_data(data_file)
    mask = data['types'] == particle_type
    positions = data['positions'][mask]
    
    if len(positions) == 0:
        print(f"No particles of type {particle_type} found.")
        return
    
    proj_pos, density_map = project_to_2d(positions, extent=extent)
    
    fig, ax = plt.subplots()
    im = ax.imshow(density_map, extent=extent, origin='lower', norm=LogNorm(),
                   cmap='viridis', interpolation='nearest')
    ax.scatter(proj_pos[:, 0], proj_pos[:, 1], s=0.5, c='red', alpha=0.5)  # Overlay particles
    ax.set_xlabel('X (kpc)')
    ax.set_ylabel('Y (kpc)')
    ax.set_title(f'2D Projection: Type {particle_type} Particles (N={len(positions):,})', color='white')
    plt.colorbar(im, ax=ax, label='Particle Density')
    
    if save:
        plt.savefig(f'{output_dir}/simulation_projection_type{particle_type}.png', dpi=300, bbox_inches='tight')
    plt.show()


def main():
    parser = argparse.ArgumentParser(description='Visualize galaxy and cluster models.')
    parser.add_argument('--model', choices=['galaxy', 'cluster', 'simulation'], default='galaxy',
                        help='Type of model to visualize')
    parser.add_argument('--input', type=str, help='Input HDF5 file for simulation')
    parser.add_argument('--radius', type=float, default=50.0, help='Model radius in kpc')
    parser.add_argument('--profile', choices=['NFW', 'beta'], default='NFW', help='Cluster profile type')
    parser.add_argument('--particle_type', type=int, default=1, help='Particle type for projection')
    parser.add_argument('--output', type=str, default='plots/', help='Output directory')
    parser.add_argument('--no-save', action='store_true', help='Do not save figures')
    
    args = parser.parse_args()
    
    if args.model == 'galaxy':
        plot_galaxy(model_type='disk', radius=args.radius, output_dir=args.output, save=not args.no_save)
    elif args.model == 'cluster':
        plot_cluster_radius(radius=args.radius, profile=args.profile, output_dir=args.output, save=not args.no_save)
    elif args.model == 'simulation':
        if not args.input:
            raise ValueError("Input file required for simulation visualization.")
        plot_simulation_projection(args.input, particle_type=args.particle_type,
                                   output_dir=args.output, save=not args.no_save)
    else:
        print("Unknown model type. Use --help for options.")


if __name__ == '__main__':
    main()