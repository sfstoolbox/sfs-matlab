# gplatex_plot_wavefield.gp
#   Usage: call 'gplatex_plot_wavefield.gp' xmin xmax ymin ymax datafile
#
#   Input parameters:
#       xmin,xmax   - xrange
#       ymin,ymax   - yrange
#       datafile    - name of the binary data file
#
#   Gnuplot version for the plot_wavefield function. the oputput will be a eps
#   and a tex file to included in any paper.
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
filename = '"$7"'
outfile = '"$8'.'.tex"'

# wxt terminal
#set terminal wxt size 350,262 enhanced font 'Verdana,10' persist
# epslatex terminal
set terminal epslatex color colortext
set output @outfile

unset key

# FIXME: the right size of the whole plot is difficulty to set, because it
# depends on the chosen X and Y values. I should check, if there exist an easy
# method to crop the produced eps file without loosing the right position of the
# text labels. Another possibility will be to calculate the splot size from the
# given X and Y values.
set size ratio -1 0.65,0.5
set bmargin screen 0.1
set tmargin screen 0.48
#set size ratio -1 0.65,0.65
#set bmargin screen 0.06
#set tmargin screen 0.64

# ---  Image plot settings ---
# Use clipping to show a reasonable range of colors
set cbrange [-1:1]
set colorbox
set palette gray

set xrange [xmin:xmax]
set yrange [ymin:ymax]
set tics scale 0.75
set border linewidth 2
set xtics 1
set ytics 1
set xlabel '$x$ (m)'
set ylabel '$y$ (m)'

# --- Drawing loudspeakers function ---
# Number of allready drawn loudspeakers
dLS = 0
# If the distance of the loudspeaker is large enough, draw the loudspeakers. If
# not, set the yrange to a larger value in order to avoid an empty area in the
# plot.
if(LSdist>0.1) call 'gp_draw_loudspeakers.gnu'; else set yrange[ymin+0.1:ymax]

plot @filename binary matrix with image
