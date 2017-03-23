#gp_set_head sets a head with nose as an gnuplot object
#   Usage: call 'gp_set_head.gnu' 'x' 'y' 'phi'
#
#   Input parameters:
#       x, y        - head position / m
#       phi         - orientation of the head / rad
#
#   gp_set_head sets two gnuplot objects in the form of a head consisting of one
#   circle and a triangle for the nose at the given position and with the nose
#   pointing in the direction of phi. for the given orientation. A phi of 0
#   means the head is pointing towards the x-axis.
#   This function counts the global variable object_number two up.
#
#   see also: gp_set_loudspeaker

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
if ($#!=3) print 'gp_set_head needs 3 input parameters'; exit

# Getting the input parameters
x = $0
y = $1
phi = $2+pi/2

# Initialize an object number
if (!exists("object_number")) object_number = 1;

# size of the head (20 cm)
s = 0.2

# Set the nose and move it to the right direction
set object object_number polygon from \
    (s-0.06)*sin(phi)+x,   -(s-0.06)*cos(phi)+y    to \
    (s-0.11)*cos(phi)+x,   (s-0.11)*sin(phi)+y     to \
    -(s-0.11)*cos(phi)+x,  -(s-0.11)*sin(phi)+y    to \
    (s-0.06)*sin(phi)+x,   -(s-0.06)*cos(phi)+y
set object object_number fc rgb 'black' fillstyle solid 0.3 lw 0.5 front
object_number = object_number+1
# Set the head
set object object_number circle at x,y size s/2.0
set object object_number fc rgb 'black' fillstyle solid 0.3 lw 0.5 front
object_number = object_number+1
