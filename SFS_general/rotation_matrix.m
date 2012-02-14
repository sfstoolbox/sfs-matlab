function R = rotation_matrix(phi,orientation)
%ROTATION_MATRIX returns a 2D rotation matrix for the given angle
%
%   Usage: R = rotation_matrix(phi)
%          R = rotation_matrix(phi,'counterclockwise')
%          R = rotation_matrix(phi,'clockwise')
%
%   Input parameters:
%       phi         - angle to rotate the given dim dimension vector (rad)
%       orientation - orientation of the rotation 
%                     (default: 'counterclockwise')
%
%   Output parameters:
%       R       - 2x2 rotation matrix to apply to your vector to rotate: 
%                 R*y
%
%   ROTATION_MATRIX(phi,orientation) returns a rotation matrix R, which is
%   able to rotate a given dim dimensional vector about phi.
%
%   see also: echo_direction
%

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
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
