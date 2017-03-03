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
# The MIT License (MIT)                                                      *
#                                                                            *
# Copyright (c) 2010-2017 SFS Toolbox Developers                             *
#                                                                            *
# Permission is hereby granted,  free of charge,  to any person  obtaining a *
# copy of this software and associated documentation files (the "Software"), *
# to deal in the Software without  restriction, including without limitation *
# the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
# and/or  sell copies of  the Software,  and to permit  persons to whom  the *
# Software is furnished to do so, subject to the following conditions:       *
#                                                                            *
# The above copyright notice and this permission notice shall be included in *
# all copies or substantial portions of the Software.                        *
#                                                                            *
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
# IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
# FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
# LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
# FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
# DEALINGS IN THE SOFTWARE.                                                  *
#                                                                            *
# The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
# methods like wave field synthesis or higher order ambisonics.              *
#                                                                            *
# http://sfstoolbox.org                                 sfstoolbox@gmail.com *
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
