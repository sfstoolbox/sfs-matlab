#!/usr/bin/gnuplot
#
#
# AUTHOR: Hagen Wierstorf
set macros
file = '"$0"'

# Set the output of the following plot to a table in order to achieve that is it
# not shown in the current terminal
set table '/dev/null'

# Function to create the right call function
add_loudspeaker(u,v,w,x) = \
    sprintf('call "gp_set_loudspeakers.gnu" "%f" "%f" "%f" "%f";',u,v,w,x);
# Initialize command string
CMD = ''
# Do a dummy plot to read the loudspeaker position data
plot @file u 1:(CMD = CMD.add_loudspeaker($$0+1,$$1,$$2,$$3))
# Execute the loudspeaker drawing command
eval(CMD)

# Restore the terminal
unset table
