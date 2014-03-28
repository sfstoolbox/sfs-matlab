function [Jsph, H2sph, Ysph, x, y, z] = sphbasis_mono_XYZgrid(X,Y,Z,f,xq,conf)
%Evaluate  basis functions for given grid in cartesian coordinates
%
%
%   Usage: [Jsph, H2sph, Ysph, x, y, z] = sphbasis_mono_XYZgrid(X,Y,Z,f,xq,conf)
%
%   Input parameters:
%       X           - x-axis / m; single value or [xmin,xmax]
%       Y           - y-axis / m; single value or [ymin,ymax]
%       Z           - z-axis / m; single value or [zmin,zmax]
%       f           - frequency in Hz
%       xq          - optional center of coordinate system
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       Jsph        - cell array of spherical bessel functions
%       H2sph       - cell array of spherical hankel functions of 2nd kind
%       Ysph        - cell array of spherical harmonics
%       x           - corresponding x axis / m
%       y           - corresponding y axis / m
%       z           - corresponding z axis / m
%
%   SPHBASIS_MONO_XYZGRID(X,Y,Z,f,xq,conf) computes spherical basis functions 
%   for given grid in cartesian coordinates. This is a wrapper function for
%   sphbasis_mono.
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
nargmin = 4;
nargmax = 6;
narginchk(nargmin,nargmax);
isargvector(X,Y,Z);
isargpositivescalar(f);
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

% Create a x-y-grid
[xx,yy,zz,x,y,z] = xyz_grid(X,Y,Z,conf);

k = 2*pi*f/conf.c;  % wavenumber

% shift coordinates to expansion coordinate
xx = xx-xq(1);
yy = yy-xq(2);
zz = zz-xq(3);

% coordinate transformation
r = vector_norm(cat(3,xx,yy,zz),3);
phi = atan2(yy,xx);
theta = asin(zz./r);

[Jsph, H2sph, Ysph] = sphbasis_mono(r, theta, phi, k, conf);

end




