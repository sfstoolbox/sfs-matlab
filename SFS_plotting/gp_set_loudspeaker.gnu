#gp_set_loudspeaker sets a loudspeaker as an gnuplot object
#   Usage: call 'gp_set_loudspeaker.gnu' 'x0' 'y0' 'phi' 'activity' 'lssize'
#
#   Input parameters:
#       x0, y0      - loudspeaker position / m
#       phi         - orientation of the loudspeaker / rad
#       activity    - activity of the loudspeaker (0..1)
#       lssize      - size of the loudspeaker / m
#
#   gp_set_loudspeaker sets a single gnuplot object in the
#   form of a loudspeaker at the given position and for the given orientation.
#   The activity gives the color of the speaker ranging from 0 (white) to 1
#   (gray). This file is run in the loop gp_draw_loudspeakers.
#   For every added loudspeaker the global variable object_number is counted one
#   up and is accessable in your gnuplot code.
#
#   see also: gp_draw_loudspeakers, gp_set_head

#*****************************************************************************
# Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
#                         Assessment of IP-based Applications                *
#                         Telekom Innovation Laboratories, TU Berlin         *
#                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
#                                                                            *
# Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
#                         Universitaet Rostock                               *
#                         Richard-Wagner-Strasse 31, 18119 Rostock           *
#                                                                            *
# This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
#                                                                            *
# The SFS is free software:  you can redistribute it and/or modify it  under *
# the terms of the  GNU  General  Public  License  as published by the  Free *
# Software Foundation, either version 3 of the License,  or (at your option) *
# any later version.                                                         *
#                                                                            *
# The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
# WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
# FOR A PARTICULAR PURPOSE.                                                  *
# See the GNU General Public License for more details.                       *
#                                                                            *
# You should  have received a copy  of the GNU General Public License  along *
# with this program.  If not, see <http://www.gnu.org/licenses/>.            *
#                                                                            *
# The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
# field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
# ambisonics.                                                                *
#                                                                            *
# http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
#*****************************************************************************


# Checking if we have enough input parameters
if ('$#'!=5) { print 'gp_set_loudspeakers needs 5 input parameters'; exit }

# Getting the input parameters
x0 = $0
y0 = $1
p = $2
activity = $3
lssize = $4

# Initialize an object number
if (!exists("object_number")) { object_number = 1; }

# Fixing loudspeaker size (because we draw a line around the loudspeaker
a = lssize-0.01;

# Set the loudspeaker at the given position and rotate it by a rotation matrix:
set object object_number polygon from \
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
set object object_number fc rgb "black" fillstyle solid 0.5*activity lw 1 front
object_number = object_number+1
