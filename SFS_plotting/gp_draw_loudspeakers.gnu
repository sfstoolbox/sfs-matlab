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
