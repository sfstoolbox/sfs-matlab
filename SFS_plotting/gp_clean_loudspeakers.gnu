# gp_clean_loudspeakers
#
# Function to unset the pbjects and arrows used to draw the loudspeakers
#
# AUTHOR: Hagen Wierstorf

if (!exists("dLS")) dLS = 0;

unset for [i=dLS+1:nLS] object i
unset for [i=dLS+1:nLS] arrow i

