function m = generate_colormap(table,n)
%GENERATE_COLORMAP creates a Matlab colormap from the given table
%
%   usage: m = generate_colormap(table,n)
%
%   Input parameters:
%       table - matrix containing the color values as columns [r g b]. The
%               single colors go from 0 to 255
%       n     - size of the colormap
%
%   Output parameters:
%       m     - colormap
%
%   GENERATE_COLORMAP(table,n) returns an n-by-3 matrix containing a colormap.
%   The color values are specified in table, which will be interpolated to the
%   desired number of entries n.
%
% See also: moreland, chromajs

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


%% ===== Checking input parameters =======================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);


%% ===== Computation =====================================================
m = zeros(n,3);
for ii=1:n
    for jj=1:3
        m(ii,jj)=interp1(linspace(1,n,size(table,1)),table(:,jj),ii);
    end
end
m=m/256;
