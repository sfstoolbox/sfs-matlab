function S = line_source(x,y,xs,f,conf)
%LINE_SOURCE returns the Green's function for a line source
%   Usage: S = line_source(x,y,xs,ys,omega,conf)
%          S = line_source(x,y,xs,ys,omega)
%
%   Input options:
%       x,y      - x,y points for which the Green's function should be calculated
%       xs,ys    - position of the point source
%       f        - frequency of the point source
%       conf     - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       S        - Wave field of a point source located at xs,ys
%
%   LINE_SOURCE(x,y,xs,f) calculates the wave field of a line source
%   located at xs,ys for the given points x,y and the frequency f. The wave
%   field is calculated by the Greens function.
%
%   References:
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: point_source

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 5;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(x,y);
xs = position_vector(xs);
isargpositivescalar(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
c = conf.c;


%% ===== Computation =====================================================
omega = 2*pi*f;
% Source model for a line source: 2D Green's function.
%
%              i   (1)/ w        \
% G(x-xs,w) =  -  H0  | - |x-xs| |
%              4      \ c        /
%
% see: Williams1999, p. 266
%
S = 1i/4 * besselh(0,1,omega/c* ...
    sqrt( (x-xs(1)).^2 + (y-xs(2)).^2 ));
