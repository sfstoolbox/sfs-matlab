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
% Made both matrices the same size
m1 = size(x1,1);
m2 = size(x2,1);
m = max(m1,m2);
n = size(x1,2);
if m1>m2
    x2 = repmat(x2,[m 1]);
elseif m2>m1
    x1 = repmat(x1,[m 1]);
end
% Calculate direction vectors
directions = zeros([m n]);
for ii=1:size(x1,1)
    directions(ii,:) = (x2(ii,:)-x1(ii,:)) / norm(x2(ii,:)-x1(ii,:));
end
