#gp_draw_loudspeakers adds loudspeaker after their positions given in a file
#
#   Usage: call 'gp_draw_loudspeakers.gnu' 'pos_file.txt' ['lssize']
#
#   Input parameters:
#       file    - file in which the positions of the loudspeaker are stored. The
#                 file has to had the following columns: x0, y0, phi, activity
#       lssize  - optional argument given the overall size of the loudspeaker.
#                 The default size is 0.1
#
#   gp_draw_loudspeakers adds a bunch of loudspeakers to your plot. The
#   number of speakers depends on the lines given in the file and places them
#   accordingly to the x0, y0 coordinates. phi gives the rotation of the speaker
#   and activity the color ranging from white (0) to gray (1).
#
#   see also: gp_set_loudspeakers

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


# getting the input parameters
set macros
file = "'$0'"
if ('$#'==2) { lssize = $1; } else { lssize = 0.1 }

# Set the output of the following plot to a table in order to achieve that is it
# not shown in the current terminal
set table '/dev/null'

# Function to create the right call function
add_loudspeaker(v,w,x,y) = \
    sprintf('call "gp_set_loudspeaker.gnu" "%f" "%f" "%f" "%f" "%f";',\
    v,w,x,y,lssize);
# Initialize command string
CMD = ''
# Do a dummy plot to read the loudspeaker position data
plot @file u 1:(CMD = CMD.add_loudspeaker($$1,$$2,$$3,$$4))
# Execute the loudspeaker drawing command
eval(CMD)

# Restore the terminal
unset table
