% simulates a monofrequent line source for a rectangular grid
% S.Spors 10.6.2005
function [S] = line_source(x,y,x_s,y_s,k_s)


[X,Y]=meshgrid(x,y);
r = sqrt( (X'-x_s).^2 + (Y'-y_s).^2);

S = besselh(0,2,k_s*r);