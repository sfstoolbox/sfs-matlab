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
