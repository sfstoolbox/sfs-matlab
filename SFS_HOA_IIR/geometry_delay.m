

nLS = 24;
R = 1;

deg = pi/180;
w = linspace(1,360,1000);
si = sin(w*deg);
co = cos(w*deg);



hold on
draw_loudspeakers(config.x0/config.R,config.n0,1);
plot(si,co,'k');
plot(1.5, 1,'-.ok','MarkerEdgeColor','k',...
                'MarkerFaceColor','k','MarkerSize',10)
l = quiver(0,-1.5,0,1.5,2,'k:');
m = quiver(-1.5,0,1.5,0,2,'k:');
n = quiver(0,0,0.7,0.71414,0,'k');
o = quiver(0.7,0.71414,0.8,0.28585,0,'k:');
p = quiver(0,0,1.5,1,0,'k');




axis equal
hold off
