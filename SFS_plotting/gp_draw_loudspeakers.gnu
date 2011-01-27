# gp_draw_loudspeakers

# AUTHOR: Hagen Wierstorf

if (!exists("dLS")) dLS = 0;

# x-Dimension of loudspeaker
b = 0.084
# draw the rectangular object of the loufspeaker
set object dLS+1 rect from x0_start-b/2+dLS*LSdist,-0.1 to \
    x0_start+b/2+dLS*LSdist,-0.045 front
set object dLS+1 rect fc rgb "black" fillstyle solid 1.0
# draw the triangular object of the loudspeaker
set arrow dLS+1 from x0_start+dLS*LSdist,0 to x0_start+dLS*LSdist,-0.05 front head
set arrow dLS+1 filled lc rgb "black"
set arrow dLS+1 size 0.075,45
dLS = dLS+1
if(dLS<nLS) reread; else dLS = 0
