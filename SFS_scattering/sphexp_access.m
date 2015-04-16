function a = sphexp_access(A, m1, n1, m2, n2)
%Access array elements of spherical expansion coefficients
%
%   Usage: a = sphexp_access(A, m1, n1, m2, n2)
%
%   Input parameters:
%       A           - 1D, 2D array of expansion coefficients
%       m1          - order of 1st dimension
%       n1          - degree of 1st dimension (optional, default = |m1|)
%       m2          - order of 2st dimension (optional, default = 0)
%       n2          - degree of 2st dimension (optional, default = |m2|)
%
%   Output parameters:
%       a           - coefficients
%
%   SPHEXP_ACCESS(A, m1, n1, m2, n2)
%
%   see also: sphexp_access

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
nargmin = 2;
nargmax = 5;
narginchk(nargmin,nargmax);
isargvector(m1);
if nargin == nargmin
  n1 = abs(m1);
else  
  isargvector(n1);
end
if nargin < 4
  m2 = 0;
else
  isargvector(m2);
end
if nargin<nargmax
  n2 = abs(m2);
else
  isargvector(n2);
end

%% ===== Computation ====================================================

a = zeros(length(m1)*length(n1),length(m2)*length(n2));

[m1, n1] = meshgrid(m1, n1);
[m2, n2] = meshgrid(m2, n2);

s1 = abs(m1(:)) <= n1(:);
s2 = abs(m2(:)) <= n2(:);

if any(s1) && any(s2)
  [l1, l2] = sphexp_index(m1(s1), n1(s1), m2(s2), n2(s2));
  a(s1,s2) = A(l1,l2);
end

end
