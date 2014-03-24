function G = greens_function_mono(x,y,z,xs,src,f,conf)
%GREENS_FUNCTION_MONO returns a Green's function in the frequency domain
%
%   Usage: G = greens_function_mono(x,y,z,xs,src,f,[conf])
%
%   Input options:
%       x,y,z   - x,y,z points for which the Green's function should be
%                 calculated / m
%       xs      - position of the source
%       src     - source model of the Green's function. Valid models are:
%                   'ps' - point source
%                   'ls' - line source
%                   'pw' - plane wave
%       f       - frequency of the source / Hz
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       G       - Green's function evaluated at the points x,y,z
%
%   GREENS_FUNCTION_MONO(x,y,z,xs,src,f) calculates the Green's function for the
%   given source model located at xs for the given points x,y and the frequency
%   f.
%
%   References:
%       H. Wierstorf (2014) - "Perceptual Assessment of Sound Field Synthesis",
%       PhD thesis, TU Berlin
%
%   see also: sound_field_mono

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
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


%% ===== Checking of input  parameters ==================================
% disabled checking for performance reasons
nargmax = 7;
if nargin<nargmax
    conf = SFS_config;
end


%% ===== Configuration ==================================================
c = conf.c;
phase = conf.phase;


%% ===== Computation =====================================================
% frequency
omega = 2*pi*f;
% calculate Green's function for the given source model
if strcmp('ps',src)
    % Source model for a point source: 3D Green's function.
    %
    %              1  e^(-i w/c |x-xs|)
    % G(x-xs,w) = --- -----------------
    %             4pi      |x-xs|
    %
    % see: Wierstorf (2014), p.21 (2.28)
    %
    G = 1/(4*pi) * exp(-1i*omega/c .* sqrt((x-xs(1)).^2+(y-xs(2)).^2+(z-xs(3)).^2)) ./ ...
            sqrt((x-xs(1)).^2+(y-xs(2)).^2+(z-xs(3)).^2);

elseif strcmp('ls',src)
    % Source model for a line source: 2D Green's function.
    %
    %              i   (2) / w        \
    % G(x-xs,w) =  -  H0  |  - |x-xs|  |
    %              4       \ c        /
    %
    % see: Wierstorf (2014), p.22 (2.32)
    %
    G = 1i/4 * besselh(0,2,omega/c* ...
        sqrt( (x-xs(1)).^2 + (y-xs(2)).^2 + (z-xs(3)).^2 ));

elseif strcmp('pw',src)
    % Source model for a plane wave:
    %
    % G(x,w) = e^(-i w/c n x)
    %
    % see: Wierstorf (2014), p.21 (2.23)
    %
    % direction of plane wave
    nxs = xs / norm(xs);
    %
    % The following code enables us to replace this two for-loops
    % for ii = 1:size(x,1)
    %     for jj = 1:size(x,2)
    %         S(ii,jj) = exp(-1i*omega/c.*nxs*[x(ii,jj) y(ii,jj) z(ii,jj)]');
    %     end
    % end
    %
    % Get a matrix in the form of
    % 1 1 1 0 0 0 0 0 0
    % 0 0 0 1 1 1 0 0 0
    % 0 0 0 0 0 0 1 1 1
    E = eye(3*size(x,1));
    E = E(1:3:end,:)+E(2:3:end,:)+E(3:3:end,:);
    % Multiply this matrix with the plane wave direction
    N = repmat(nxs,size(x,1)) .* E;
    % Interlace x,y,z into one matrix
    % x11 x12 ... x1m
    % y11 y12 ... y1m
    % z11 z12 ... z1m
    % .   .       .
    % .   .       .
    % xn1 xn2 ... xnm
    % yn1 yn2 ... ynm
    % zn1 zn2 ... znm
    XYZ = zeros(3*size(x,1),size(x,2));
    XYZ(1:3:end,:) = x;
    XYZ(2:3:end,:) = y;
    XYZ(3:3:end,:) = z;
    % calculate sound field
    G = exp(-1i*omega/c.*N*XYZ);

else
    error('%s: %s is not a valid source model for the Green''s function', ...
        upper(mfilename),src);
end

% Add phase to be able to simulate different time steps
G = G .* exp(-1i*phase);
