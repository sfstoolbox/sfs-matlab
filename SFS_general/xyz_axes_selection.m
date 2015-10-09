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
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
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
