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
%   See also: sin, cos

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
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
% Fill up matrix to rotate around the given axis
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
