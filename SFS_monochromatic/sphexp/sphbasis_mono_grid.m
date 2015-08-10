function [jn, h2n, Ynm, x, y, z] = sphbasis_mono_grid(X,Y,Z,Nse,f,xq,conf)
%SPHBASIS_MONO_GRID(X,Y,Z,f,xq,conf) evaluates spherical basis functions for 
%given grid in cartesian coordinates
%
%   Usage: [jn, h2n, Ynm, x, y, z] = sphbasis_mono_grid(X,Y,Z,Nse,f,xq,conf)
%
%   Input parameters:
%       X           - x-axis / m; single value or [xmin,xmax] or nD-array
%       Y           - y-axis / m; single value or [ymin,ymax] or nD-array
%       Z           - z-axis / m; single value or [zmin,zmax] or nD-array
%       Nse         - maximum order of spherical basis functions
%       f           - frequency in Hz
%       xq          - optional center of coordinate system
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       jn          - cell array of spherical bessel functions
%       h2n         - cell array of spherical hankel functions of 2nd kind
%       Ynm         - cell array of spherical harmonics
%       x           - corresponding x axis / m
%       y           - corresponding y axis / m
%       z           - corresponding z axis / m
%
%   SPHBASIS_MONO_GRID(X,Y,Z,Nse,f,xq,conf) computes spherical basis functions
%   for given grid in cartesian coordinates. This is a wrapper function for
%   sphbasis_mono.
%   For the input of X,Y,Z (DIM as a wildcard) :
%     * if DIM is given as single value, the respective dimension is
%     squeezed, so that dimensionality of the simulated sound field P is
%     decreased by one.
%     * if DIM is given as [dimmin, dimmax], a linear grid for the
%     respective dimension with a resolution defined in conf.resolution is
%     established
%     * if DIM is given as n-dimensional array, the other dimensions have
%     to be given as n-dimensional arrays of the same size or as a single value.
%     Each triple of X,Y,Z is interpreted as an evaluation point in an
%     customized grid. For this option, plotting and normalisation is disabled.
%
%   see also: sphbasis_mono

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
nargmin = 5;
nargmax = 7;
narginchk(nargmin,nargmax);

isargnumeric(X,Y,Z);
% unique index encoding which dimension is an nd-array
customGrid = (numel(X) > 2) + 2*(numel(Y) > 2) + 4*(numel(Z) > 2);
switch customGrid
  case 1
    isargscalar(Y,Z);
  case 2
    isargscalar(X,Z);
  case 3
    isargequalsize(X,Y); isargscalar(Z);
  case 4
    isargscalar(X,Y);
  case 5
    isargequalsize(X,Z); isargscalar(Y);
  case 6
    isargequalsize(Y,Z); isargscalar(X);
  case 7
    isargequalsize(X,Y,Z);
  otherwise
    isargvector(X,Y,Z);
end

isargpositivescalar(f,Nse);
if nargin<nargmax
  conf = SFS_config;
else
  isargstruct(conf);
end
if nargin == nargmin
  xq = [0, 0, 0];
end
isargposition(xq);

%% ===== Computation ====================================================
if customGrid
  % just copy everything
  xx = X; yy = Y; zz = Z;
  x = X;   y = Y;  z = Z;
else
  % Create a x-y-grid
  [xx, yy, zz, x, y, z] = xyz_grid(X, Y, Z, conf);
end

k = 2*pi*f/conf.c;  % wavenumber

% shift coordinates to expansion coordinate
xx = xx - xq(1);
yy = yy - xq(2);
zz = zz - xq(3);

% coordinate transformation
r = sqrt(xx.^2 + yy.^2 + zz.^2);
phi = atan2(yy, xx);
theta = asin(zz./r);

[jn, h2n, Ynm] = sphbasis_mono(r, theta, phi, Nse, k, conf);
