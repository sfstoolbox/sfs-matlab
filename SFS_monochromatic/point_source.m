function S = point_source(x,y,xs,ys,f,conf)
%POINT_SOURCE returns the Green's function for a point source
%   Usage: S = point_source(x,y,x0,y0,omega,conf)
%          S = point_source(x,y,x0,y0,omega)
%
%   Input options:
%       x,y     - x,y points for which the Green's function should be calculated
%       xs,ys   - position of the point source
%       f       - frequency of the point source
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       S       - Wave field of a point source located at x0,y0
%
%   POINT_SOURCE(x,y,xs,ys,f) calculates the wave field of a point source
%   located at xs,ys for the given points x,y and the frequency f. The wave
%   field is calculated by the Greens function.
%
%   References:
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: line_source

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 6;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix({x,y},{'x','y'});
isargscalar({xs,ys},{'xs','ys'});
isargpositivescalar({f},{'f'});
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct({conf},{'conf'});
end


%% ===== Configuration ==================================================
c = conf.c;


%% ===== Computation =====================================================
omega = 2*pi*f;
% Source model for a point source: 3D Green's function.
%
%              1  e^(-i w/c |x-xs|)
% G(x-xs,w) = --- -----------------
%             4pi      |x-xs|
%
% see: Williams1999, p. 198
%
S = 1/(4*pi) * exp(1i*omega/c.*sqrt((x-xs).^2+(y-ys).^2)) ./ ...
        sqrt((x-xs).^2+(y-ys).^2);
