function [jn, h2n, Ynm, x, y, z] = sphbasis_mono_grid(X,Y,Z,Nse,f,xq,conf)
%SPHBASIS_MONO_GRID(X,Y,Z,f,xq,conf) evaluates spherical basis functions for 
%given grid in cartesian coordinates
%
%   Usage: [jn, h2n, Ynm] = sphbasis_mono_grid(X,Y,Z,Nse,f,xq,conf)
%
%   Input parameters:
%       X           - x-axis / m; single value or [xmin,xmax] or nD-array
%       Y           - y-axis / m; single value or [ymin,ymax] or nD-array
%       Z           - z-axis / m; single value or [zmin,zmax] or nD-array
%       Nse         - maximum order of spherical basis functions
%       f           - frequency in Hz
%       xq          - center of coordinate system
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       jn          - cell array of spherical bessel functions
%       h2n         - cell array of spherical hankel functions of 2nd kind
%       Ynm         - cell array of spherical harmonics
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
% Copyright (c) 2010-2016 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2016 Institut fuer Nachrichtentechnik                   *
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
nargmin = 7;
nargmax = 7;
narginchk(nargmin,nargmax);

isargnumeric(X,Y,Z);
isargpositivescalar(Nse,f);
isargcoord(xq);

%% ===== Computation ====================================================
[x,y,z] = xyz_grid(X, Y, Z, conf);

k = 2*pi*f/conf.c;  % wavenumber

% shift coordinates to expansion coordinate
x = x - xq(1);
y = y - xq(2);
z = z - xq(3);
% coordinate transformation
r = sqrt(x.^2 + y.^2 + z.^2);
phi = atan2(y, x);
theta = asin(z./r);

[jn, h2n, Ynm] = sphbasis_mono(r, theta, phi, Nse, k, conf);