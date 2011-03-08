#gp_set_loudspeakers sets a loudspeaker as an gnuplot object
#   Usage: call 'gp_set_loudspeakers.gnu' 'n' 'x0' 'y0' 'phi' 'activity' 'lssize'
#
#   Input parameters:
#       n           - number of object
#       x0, y0      - loudspeaker position
#       phi         - orientation of the loudspeaker
#       activity    - activity of the loudspeaker (0..1)
#       lssize      - size of the loudspeaker
#
#   gp_set_loudspeakers sets a single gnuplot object with the number n in the
#   form of a loudspeaker at the given position and for the given orientation.
#   The activity gives the color of the speaker ranging from 0 (white) to 1
#   (gray). This file is ruin in the loop gp_draw_loudspeakers
#
#   see also: gp_draw_loudspeakers

# AUTHOR: Hagen Wierstorf

# Checking if we have enough input parameters
if ($#!=6) print 'gp_set_loudspeakers needs 6 six input parameters'; exit

# Getting the input parameters
n = $0
x0 = $1
y0 = $2
p = $3+pi/2.0
activity = $4
lssize = $5

# Fixing loudspeaker size (because we draw a line around the loudspeaker
a = lssize-0.01;

# Set the loudspeaker at the given position and rotate it by a rotation matrix:
set object n polygon from \
-a*cos(p)+a/2*sin(p)+x0,   -a*sin(p)-a/2*cos(p)+y0    to \
-a*cos(p)-a/2*sin(p)+x0,   -a*sin(p)+a/2*cos(p)+y0    to \
-a/2*cos(p)-a/2*sin(p)+x0, -a/2*sin(p)+a/2*cos(p)+y0  to \
-a/2*cos(p)-a/6*sin(p)+x0, -a/2*sin(p)+a/6*cos(p)+y0  to \
0*cos(p)-a/2*sin(p)+x0,    0*sin(p)+a/2*cos(p)+y0     to \
0*cos(p)+a/2*sin(p)+x0,    0*sin(p)-a/2*cos(p)+y0     to \
-a/2*cos(p)+a/6*sin(p)+x0, -a/2*sin(p)-a/6*cos(p)+y0  to \
-a/2*cos(p)+a/2*sin(p)+x0, -a/2*sin(p)-a/2*cos(p)+y0  to \
-a*cos(p)+a/2*sin(p)+x0,   -a*sin(p)-a/2*cos(p)+y0
# Set the color etc.
set object n fc rgb "black" fillstyle solid 0.5*activity lw 1 front
