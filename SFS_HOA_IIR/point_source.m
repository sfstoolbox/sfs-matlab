% simulates a monofrequent 3D point source for a rectangular grid
% S.Spors 24.3.2004
function [S] = point_source(x,y,x_s,y_s,k_s)


[X,Y]=meshgrid(x,y);
r = sqrt( (X'-x_s).^2 + (Y'-y_s).^2);

S = 1./r .* exp(-1i*k_s*r);