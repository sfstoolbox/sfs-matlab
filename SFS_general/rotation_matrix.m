function R = rotation_matrix(phi,orientation)
%ROTATION_MATRIX returns a 2D rotation matrix for the given angle
%
%   Usage: R = rotation_matrix(phi,[orientation])
%
%   Input parameters:
%       phi         - angle to rotate the given dim dimension vector (rad)
%       orientation - orientation of the rotation, 'clockwise' or
%                     'counterclockwise' (default: 'counterclockwise')
%
%   Output parameters:
%       R       - 2x2 rotation matrix to apply to your vector to rotate:
%                 R*y
%
%   ROTATION_MATRIX(phi,orientation) returns a rotation matrix R, which is
%   able to rotate a given dim dimensional vector about phi.
%
%   see also: echo_direction

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
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
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
isargscalar(phi)
if nargin < 2
    % Set defualt orientation of the rotation
    orientation = 'counterclockwise';
else
    isargchar(orientation);
end


%% ===== Computation ====================================================
% Rotation matrix (see: http://en.wikipedia.org/wiki/Rotation_matrix)
switch orientation
    case 'counterclockwise'
        R = [cos(phi) -sin(phi); ...
             sin(phi) cos(phi)];
    case 'clockwise'
        R = [cos(phi) sin(phi); ...
             -sin(phi) cos(phi)];
end
