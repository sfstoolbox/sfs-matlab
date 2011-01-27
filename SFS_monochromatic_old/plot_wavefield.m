% plot reproduced wave field
% S.Spors / 27.7.2007

P = 0.5*P/abs(P(end/2,end/2));
imagesc(x,y,real(P'));
%imagesc(x,y,db(abs(P')));
colorbar;
caxis([-1 1]);
draw_loudspeakers(x0,n0,sign(abs(D)));
%draw_loudspeakers(x0,n0,abs(D)/max(abs(D(:))));
%axis square;
axis([x(1) x(end) y(1) y(end)]);
%axis equal;
%axis tight;
turn_imagesc;

xlabel('x -> [m]');
ylabel('y -> [m]');