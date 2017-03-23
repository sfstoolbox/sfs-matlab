function directions = direction_vector(x1,x2)
%DIRECTION_VECTOR return unit vector(s) pointing from x1 to x2
%
%   Usage: n = direction_vector(x1,x2)
%
%   Input parameters:
%       x1  - starting point(s) [1xn] or [mxn]
%       x2  - ending point(s)   [1xn] or [mxn]
%
%   Output parameters:
%       n   - unit vector(s) pointing in the direction(s) from x1 to x2
%
%   DIRECTION_VECTOR(x1,x2) calculates the unit vectors pointing from
%   n-dimensional points x1 to the n-dimensional points x2.
%
%   See also: secondary_source_positions

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


%% ===== Checking of input  parameters ===================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
if size(x1,2)~=size(x2,2)
    error('%s: x1 and x2 need to have the same dimension.',upper(mfilename));
end
if size(x1,1)~=size(x2,1) && ~(size(x1,1)==1 | size(x2,1)==1)
    error(['%s: x1 and x2 need to have the same size, or one needs to ', ...
           'be a vector.'],upper(mfilename));
end


%% ==== Main =============================================================
% Calculate direction vectors
directions = bsxfun(@minus, x2, x1);  % automatically adjusts matrix sizes
directions = bsxfun(@rdivide, directions, sqrt(sum(abs(directions).^2,2)));
