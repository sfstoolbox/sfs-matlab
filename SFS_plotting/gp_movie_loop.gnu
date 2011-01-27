# Plots the wave field movie frames (files) in an endless loop 
#
# AUTHOR: Hagen Wierstorf

set output tmpdir.'/'.frame.'.png'
plot tmpdir.'/'.frame.'.dat' binary matrix with image
frame = frame+1
if(frame<35) reread

# Endless loop for wxt terminal
#if(frame>35) frame = 12
#reread


