function [dimensions,x1,x2,x3] = xyz_axes_selection(x,y,z)
%XYZ_AXES_SELECTION returns the first non-singleton axes and a vector
%indicating which axes are selected
%
%   Usage: [dimensions,x1,x2,x3] = xyz_axes_selection(x,y,z)
%
%   Input options:
%       x,y,z      - vectors/matrices containing the x-, y- and z-axis values / m
%
%   Output options:
%       dimensions - 1x3 vector containing 1 or 0 to indicate the activity
%                    of the single dimensions in the order [x y z]
%       x1         - vector/matrix containing the first axis / m
%       x2         - vector/matrix containing the second axis / m
%       x3         - vector/matrix containing the third axis / m
%
%   XYZ_AXES_SELECTION(x,y,z) returns a indicating vector for the x-, y- and
%   z-axis if we have any activity on this axis or if it is a singleton axis.
%   In addition, the axes are reordered starting first with the non-singleton
%   axes.
%
%   See also: plot_sound_field, xyz_grid, is_dim_custom

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


%% ===== Checking of input parameters ====================================
nargmin = 3;
nargmax = 3;
narginchk(nargmin,nargmax);
isargnumeric(x,y,z);


%% ===== Computation =====================================================
dims = {x,y,z};
dimensions = ~is_dim_singleton(dims{:});
Nd = sum(dimensions);

newdims = {x(1),y(1),z(1)};  % default case, if all dimensions are singleton
newdims(Nd+1:end) = newdims(~dimensions); % move singleton dimensions to the end
newdims(1:Nd) = dims(dimensions);  % move non-singleton dimensions to the front

[x1, x2, x3] = newdims{:};
