# plot_comparison.gnu
#
# Gnuplot script to compare density profiles of galaxy/cluster models.
# Plots rho(r) vs. r (log-log scale) for NFW (dark matter halo), Hernquist (spheroid),
# and Plummer (cored stellar system) from CSV outputs.
#
# Assumes CSV format: r_kpc,rho_Msun_kpc3,M_enc_Msun,Phi_kms2
# (Columns: 1=r, 2=rho)
#
# Usage: gnuplot plot_comparison.gnu
#
# Outputs: comparison_profiles.png and comparison_profiles.pdf
#
# References:
# - NFW: Navarro et al. (1996)
# - Hernquist: Hernquist (1990)
# - Plummer: Plummer (1911)
#

# Set terminal and output
set terminal pngcairo enhanced font 'Arial,12' linewidth=2 size 1200,800
set output 'comparison_profiles.png'

# Alternative for PDF
# set terminal pdfcairo enhanced color font 'Arial,12'
# set output 'comparison_profiles.pdf'

# General settings
set logscale x 2 y  # Log-log for profiles (r from 0.1-100 kpc, rho dynamic)
set xlabel 'Radius (kpc)'
set ylabel 'Density (M_⊙ kpc^{-3})'
set title 'Comparison of Density Profiles: NFW, Hernquist, Plummer'
set key top right box opaque
set grid minor lc rgb '#dddddd'
set border 1 lw 2

# X/Y ranges (adjust based on data)
set xrange [0.1:100]
set yrange [1e6:1e12]  # Typical for 1e12-1e15 Msun systems

# Plot density (column 1=r, 2=rho) from CSV files
# Assume files in ../examples/ or current dir; adjust paths as needed
plot \
  'nfw_profile.csv' using 1:2 title 'NFW (DM Halo)' with lines lw 3 lc rgb 'blue' dashtype 1, \
  'hernquist_profile.csv' using 1:2 title 'Hernquist (Spheroid)' with lines lw 3 lc rgb 'red' dashtype 2, \
  'plummer_profile.csv' using 1:2 title 'Plummer (Cored)' with lines lw 3 lc rgb 'green' dashtype 3

# Reset output for interactive viewing
set terminal wxt
set output

# Optional: Plot enclosed mass comparison (uncomment for M(<r))
# set yrange [1e9:1e15]
# set ylabel 'Enclosed Mass (M_⊙)'
# plot \
#   'nfw_profile.csv' using 1:3 title 'NFW M(<r)' with lines lw 2 lc rgb 'blue', \
#   'hernquist_profile.csv' using 1:3 title 'Hernquist M(<r)' with lines lw 2 lc rgb 'red', \
#   'plummer_profile.csv' using 1:3 title 'Plummer M(<r)' with lines lw 2 lc rgb 'green'

# End of script