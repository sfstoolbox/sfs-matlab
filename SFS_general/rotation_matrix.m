function R = rotation_matrix(phi,dim,orientation)
%ROTATION_MATRIX returns a 3D rotation matrix for the given angle and dimension
%
%   Usage: R = rotation_matrix(phi,[dim,[orientation]])
%
%   Input parameters:
%       phi         - angle to rotate the given dim dimension vector / rad
%       dim         - dimension to turn around (default: 3, z-axis)
%       orientation - orientation of the rotation, 'clockwise' or
%                     'counterclockwise' (default: 'counterclockwise')
%
%   Output parameters:
%       R       - 3x3 rotation matrix to apply to your vector to 
%                 rotate: R*y
%                 
%
%   ROTATION_MATRIX(phi,dimension,orientation) returns a rotation matrix R, 
%   which is able to rotate a vector around the given dimension about phi.
%
%   see also: echo_direction

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
nargmin = 1;
nargmax = 3;
narginchk(nargmin,nargmax);
isargscalar(phi)
if nargin<3
    % Set defualt orientation of the rotation
    orientation = 'counterclockwise';
end
if nargin<2
    % set default rotation dimension to z-axis
    dim = 3;
end
isargchar(orientation);
isargpositivescalar(dim);


%% ===== Computation ====================================================
% Rotation matrix (see: http://en.wikipedia.org/wiki/Rotation_matrix)
% get single matrix entries
r1 = cos(phi);
r4 = cos(phi);
if strcmp('counterclockwise',orientation)
    r2 = -sin(phi);
    r3 =  sin(phi);
elseif strcmp('clockwise',orientation)
    r2 =  sin(phi);
    r3 = -sin(phi);
else
    error('%s: the given orientation "%s" is not known.', ...
        upper(mfilename),orientation);
end
% fill up matrix to rotate around the given axis
if dim==1

    R = [1 0  0;  ...
         0 r1 r2; ...
         0 r3 r4];
elseif dim==2
    R = [r1 0 r2; ...
         0  1 0;  ...
         r3 0 r4];
elseif dim==3
    R = [r1 r2 0; ...
         r3 r4 0; ...
         0  0  1];
else
    error('%s: dim has to be 1,2, or 3 and not %i',upper(mfilename),dim);
end
