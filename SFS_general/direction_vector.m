function n = direction_vector(x1,x2)
%DIRECTION_VECTOR returns a unit vector pointing from x1 to x2
%
%   Usage: n = direction_vector(x1,x2)
%
%   Input parameters:
%       x1  - starting point
%       x2  - ending point
%
%   Output parameters:
%       n   - unit vector pointing in the direction from x1 to x2
%
%   DIRECTION_VECTOR(x1,x2) calculates the unit vector pointing from the
%   n-dimensional point x1 to the n-dimensional point x2. The vectors x1
%   and x2 can each be stored in a matrix containing m different direction
%   vectors.
%
%   see also: secondary_source_positions

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


%% ===== Checking of input  parameters ===================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
if size(x1)~=size(x2)
    error('%s: x1 and x2 had to have the same size.',upper(mfilename));
end


%% ==== Main =============================================================
n = zeros(size(x1));
for ii=1:size(x1,1)
    n(ii,:) = (x2(ii,:)-x1(ii,:)) / norm(x2(ii,:)-x1(ii,:));
end
