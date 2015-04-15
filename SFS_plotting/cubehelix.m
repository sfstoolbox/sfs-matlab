function map = cubehelix(nlev,start,rots,hue,gamma)
%CUBEHELIX returns the cubehelix colormap
%
%   Usage: cubehelix([[n],start,rots,hue,gamma])
%
%   Input parameters:
%       n     - optional length of colormap (default uses the figure default
%               length)
%       start - colour to begin at (1=red, 2=green, 3=red; e.g. 0.5=purple),
%               default: 0.5
%       rots  - number of rotations, default: -1.5
%       hue   - hue intensity scaling, default: 1.2
%       gamma - intensity correction, default: 1.0
%
%   Output parameters:
%       map   - colormap [n 3]
%
%   CUBEHELIX(nlev,start,rots,hue,gamma) calculates a "cube helix" color map.
%   The colors are a tapered helix around the diagonal of the RGB colour cube,
%   from black [0,0,0] to white [1,1,1]. Deviations away from the diagonal vary
%   quadratically, increasing from zero at black to a maximum and then
%   decreasing to zero at white, all the time rotating in color.
%   The routine returns an n-by-3 matrix.
%
%   For example, to reset the colormap of the current figure:
%       colormap(cubehelix)
%
%   See arXiv:1108.5083 and https://www.mrao.cam.ac.uk/~dag/CUBEHELIX/ for more
%   details.
%
%   Original written by Dave Green and Philip Graff
%
%   See also: chromajs, moreland, plot_sound_field

%*****************************************************************************
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking input parameters =======================================
nargmin = 0;
nargmax = 5;
narginchk(nargmin,nargmax);
if nargin==0
    n = size(get(gcf,'colormap'),1);
end
if nargin==0 | nargin==1
    start = 0.5;
    rots = -1.5;
    hue = 1.2;
    gamma = 1.0;
else
    error('%s: you have to use 0, 1, or 5 input args.',upper(mfilename));
end
map = zeros(n,3);
A = [-0.14861,1.78277;-0.29227,-0.90649;1.97294,0];
for ii=1:n
    fract = (ii-1)/(n-1);
    angle = 2*pi*(start/3+1+rots*fract);
    fract = fract^gamma;
    amp = hue*fract*(1-fract)/2;
    map(ii,:) = fract+amp*(A*[cos(angle);sin(angle)])';
    for jj=1:3
        if map(ii,jj)<0
            map(ii,jj) = 0;
        elseif map(ii,jj)>1
            map(ii,jj) = 1;
        end
    end
end
