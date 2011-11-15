#gp_set_head sets a head with nose as an gnuplot object
#   Usage: call 'gp_set_head.gnu' 'x' 'y' 'phi'
#
#   Input parameters:
#       x, y        - head position
#       phi         - orientation of the head (rad)
#
#   gp_set_head sets two gnuplot objects in the form of a head consisting of one
#   circle and a triangle for the nose at the given position and with the nose
#   pointing in the direction of phi. for the given orientation. A phi of 0
#   means the head is pointing towards the x-axis.
#   This function counts the global variable object_number two up.
#
#   see also: gp_set_loudspeaker

# AUTHOR: Hagen Wierstorf

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
