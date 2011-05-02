% plot reproduced wave field
% S.Spors / 27.7.2007


P = 0.7*P./abs(P(end/2,end/2));
imagesc(x,y,real(P'));
%imagesc(x,y,db(abs(P')));
colorbar;
caxis([-1 1]);

draw_loudspeakers(x0,n0,sign(abs(D)));
turn_imagesc;

xlabel('x -> [m]');
ylabel('y -> [m]');