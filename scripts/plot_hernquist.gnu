# plot_hernquist.gnu
#
# Gnuplot script to visualize the Hernquist density profile (Hernquist 1990),
# used for modeling spheroidal galaxy components in clusters.
# Plots rho(r) and M(<r) vs. r (log-log scale) from CSV output.
#
# Assumes CSV format: r_kpc,rho_Msun_kpc3,M_enc_Msun,Phi_kms2
# (Columns: 1=r, 2=rho, 3=M_enc)
#
# Usage: gnuplot plot_hernquist.gnu
#
# Outputs: hernquist_profiles.png and hernquist_profiles.pdf
#
# References:
# Hernquist, L. 1990, ApJ, 356, 359
#

# Set terminals and outputs
set terminal pngcairo enhanced font 'Arial,12' linewidth=2 size 1200,600
set output 'hernquist_profiles.png'

# Alternative for PDF
# set terminal pdfcairo enhanced color font 'Arial,12'
# set output 'hernquist_profiles.pdf'

# Multi-panel layout: Density and Enclosed Mass
set multiplot layout 1,2

# Panel 1: Density rho(r) vs. r (log-log)
set logscale x y
set xlabel 'Radius (kpc)'
set ylabel 'Density (M_⊙ kpc^{-3})'
set title 'Hernquist Density Profile ρ(r)'
set key off
set grid minor lc rgb '#dddddd'
set border 1 lw 2
set xrange [0.1:100]
set yrange [1e6:1e12]  # Adjust based on M_tot=1e12 Msun, a=5 kpc

plot 'hernquist_profile.csv' using 1:2 with lines lw 3 lc rgb 'blue' title 'ρ(r)'

# Panel 2: Enclosed Mass M(<r) vs. r (log-log)
set logscale y
set logscale x
set xlabel 'Radius (kpc)'
set ylabel 'Enclosed Mass (M_⊙)'
set title 'Hernquist Enclosed Mass M(<r)'
set key top left box opaque
set grid minor lc rgb '#dddddd'
set border 1 lw 2
set xrange [0.1:100]
set yrange [1e9:1e12]

plot 'hernquist_profile.csv' using 1:3 with lines lw 3 lc rgb 'red' title 'M(<r)'

# Optional: Add potential Phi(r) plot (uncomment and adjust multiplot to 1,3)
# set logscale x
# set nologscale y
# set yrange [-1e5:0]  # Typical km^2/s^2
# set ylabel 'Potential Φ(r) (km² s^{-2})'
# set title 'Hernquist Potential Φ(r)'
# plot 'hernquist_profile.csv' using 1:4 with lines lw 3 lc rgb 'green' title 'Φ(r)'

# Reset for interactive
unset multiplot
set terminal wxt
set output

# End of script