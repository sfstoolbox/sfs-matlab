#!/usr/bin/gnuplot
#
#
# AUTHOR: Hagen Wierstorf
set macros
file = '"$0"'

# Store current terminal and choose a dummy terminal
set term push
set term dumb
#set output '/dev/null'

# Function to create the right call function
f(u,v,w,x) = \
    sprintf('call "gp_set_loudspeakers.gnu" "%f" "%f" "%f" "%f";',u,v,w,x);
# Initialize command string
CMD = ''
# Do a dummy plot to read the loudspeaker position data
plot @file u 1:(CMD = CMD.f($$0+1,$$1,$$2,$$3))
# Execute the loudspeaker drawing command
eval(CMD)

# Restore the terminal
set term pop
