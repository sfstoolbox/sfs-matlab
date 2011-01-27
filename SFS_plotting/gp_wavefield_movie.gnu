# gp_plot_wavefield.gp
#   Usage: call 'gp_plot_wavefield.gp' xmin xmax ymin ymax datafile
#
#   Input parameters:
#       xmin,xmax   - xrange
#       ymin,ymax   - yrange
#       L	    - array length
#       LSdist	    - distance between loudspeakers
#       X0          - center of loudspeaker array
#       datafile    - name of the binary data file
#
#   Gnuplot version for the plot_wavefield function.
#   This file has to be called by the call command!

# AUTHOR: Hagen Wierstorf

reset
set macros

# Parse input parameter
xmin = $0
xmax = $1
ymin = $2
ymax = $3
nLS = $4
LSdist = $5
x0_start = $6
#filename = '$7'
tmpdir = '$7'

# wxt terminal
#set terminal wxt size 700,524 enhanced font 'Verdana,14' persist
# png terminal
set terminal pngcairo size 700,524 enhanced font 'Verdana,14'

unset key

set size ratio -1
set border linewidth 2

# ---  Image plot settings ---
# Use clipping to show a reasonable range of colors
set cbrange [-1:1]
set colorbox
set palette gray

set xrange [xmin:xmax]
set yrange [ymin:ymax]
set tics scale 0.75
set xtics 1
set ytics 1
set xlabel 'x (m)'
set ylabel 'y (m)'

# --- Drawing loudspeakers function ---
# Number of drawn loudspeakers
dLS = 0
if(LSdist>0.1) call 'gp_draw_loudspeakers.gnu'

# Start with frame 11 in order to avoid sort problem with numbers 1-9
frame = 11
# --- 
call 'gp_movie_loop.gnu'
