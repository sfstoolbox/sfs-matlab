#gp_draw_loudspeakers adds loudspeaker after their positions given in a file
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

# AUTHOR: Hagen Wierstorf

# getting the input parameters
set macros
file = '"$0"'
if ($#==2) lssize = $1; else lssize = 0.1

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
