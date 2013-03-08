function S = plane_wave(x,y,z,xs,f,conf)
%PLANE_WAVE returns the source model for a plane wave
%
%   Usage: S = plane_wave(x,y,z,x0,omega,[conf])
%
%   Input options:
%       x,y,z   - x,y,z points for which the Green's function should be calculated
%       xs      - direction of the plane wave
%       f       - frequency of the point source
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       S       - Wave field of a plane wave traveling in the direction xs
%
%   PLANE_WAVE(x,y,z,xs,f) calculates the wave field of a plane wave going into
%   the direction xs for the given points x,y,z and the frequency f.
%
%   References:
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: point_source, line_source

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 6;
narginchk(nargmin,nargmax);
isargmatrix(x,y,z);
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
% Source model for a plane wave:
%
% S(x,w) = e^(-i w/c n x)
%
% see: Williams1999, p. 21
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
S = exp(-1i*omega/c.*N*XYZ);
