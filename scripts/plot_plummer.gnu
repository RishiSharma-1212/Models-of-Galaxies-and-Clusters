# plot_plummer.gnu
#
# Gnuplot script to visualize the Plummer density profile (Plummer 1911),
# used for modeling cored stellar systems like globular clusters in galaxy simulations.
# Plots rho(r), M(<r), and Phi(r) vs. r (appropriate scales) from CSV output.
#
# Assumes CSV format: r_kpc,rho_Msun_kpc3,M_enc_Msun,Phi_kms2
# (Columns: 1=r, 2=rho, 3=M_enc, 4=Phi)
#
# Usage: gnuplot plot_plummer.gnu
#
# Outputs: plummer_profiles.png and plummer_profiles.pdf
#
# References:
# Plummer, H. C. 1911, MNRAS, 71, 460
#

# Set terminals and outputs
set terminal pngcairo enhanced font 'Arial,12' linewidth=2 size 1200,800
set output 'plummer_profiles.png'

# Alternative for PDF (uncomment for vector output)
# set terminal pdfcairo enhanced color font 'Arial,12'
# set output 'plummer_profiles.pdf'

# Multi-panel layout: 1 row, 3 columns for rho, M(<r), Phi
set multiplot layout 1,3

# Panel 1: Density rho(r) vs. r (log-log, flat core)
set logscale x y
set xlabel 'Radius (kpc)'
set ylabel 'Density (M_⊙ kpc^{-3})'
set title 'Plummer Density Profile ρ(r)'
set key off
set grid minor lc rgb '#dddddd'
set border 1 lw 2
set xrange [0.01:20]  # Adjusted for a=1 kpc
set yrange [1e5:1e10]  # For M_tot=1e10 Msun

plot 'plummer_profile.csv' using 1:2 with lines lw 3 lc rgb 'blue' title 'ρ(r)'

# Panel 2: Enclosed Mass M(<r) vs. r (log-log)
set logscale x y
set xlabel 'Radius (kpc)'
set ylabel 'Enclosed Mass (M_⊙)'
set title 'Plummer Enclosed Mass M(<r)'
set key top left box opaque
set grid minor lc rgb '#dddddd'
set border 1 lw 2
set xrange [0.01:20]
set yrange [1e8:1e10]

plot 'plummer_profile.csv' using 1:3 with lines lw 3 lc rgb 'red' title 'M(<r)'

# Panel 3: Potential Phi(r) vs. r (linear y for depth)
set logscale x
set nologscale y
set xlabel 'Radius (kpc)'
set ylabel 'Potential Φ(r) (km² s^{-2})'
set title 'Plummer Potential Φ(r)'
set key bottom right box opaque
set grid minor lc rgb '#dddddd'
set border 1 lw 2
set xrange [0.01:20]
set yrange [-1e4:0]  # Typical for small systems

plot 'plummer_profile.csv' using 1:4 with lines lw 3 lc rgb 'green' title 'Φ(r)'

# Reset and interactive viewing
unset multiplot
set terminal wxt
set output

# End of script