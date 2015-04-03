function map = cubehelix(nlev,start,rots,hue,gamma)
% Usage:
%         colormap(CubeHelix(...))
%
%==========================================================
% Calculates a "cube helix" colour map for MATLAB. The
% colours are a tapered helix around the diagonal of the
% RGB colour cube, from black [0,0,0] to white [1,1,1].
% Deviations away from the diagonal vary quadratically,
% increasing from zero at black to a maximum and then
% decreasing to zero at white, all the time rotating in
% colour.
%
% The input values are:
%   nlev  = number of colour steps
%   start = colour to begin at (1=red, 2=green, 3=red;
%           e.g. 0.5=purple)
%   rots  = number of rotations
%   hue   = hue intensity scaling, 0=B&W
%   gamma = intensity correction
%
% The routine returns an nlev-by-3 matrix that can be used
% as a colourmap for MATLAB (function 'colormap').
%
% Use (256,0.5,-1.5,1.2,1.0) as defaults.
%
% See arXiv:1108.5083 for more details.
%
% See: https://www.mrao.cam.ac.uk/~dag/CUBEHELIX/
% This should be tested for the dB plots
%
%----------------------------------------------------------
% Original written in Fortran77 by Dave Green, 2011 Jan 10
% Transcribed to MATLAB by Philip Graff, 2011 Sept 1
%==========================================================
map=zeros(nlev,3);
A=[-0.14861,1.78277;-0.29227,-0.90649;1.97294,0];
for i=1:nlev
    fract=(i-1)/(nlev-1);
    angle=2*pi*(start/3+1+rots*fract);
    fract=fract^gamma;
    amp=hue*fract*(1-fract)/2;
    map(i,:)=fract+amp*(A*[cos(angle);sin(angle)])';
    for j=1:3
        if map(i,j)<0
            map(i,j)=0;
        elseif map(i,j)>1
            map(i,j)=1;
        end
    end
end
