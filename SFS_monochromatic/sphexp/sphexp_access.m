function a = sphexp_access(Anm,m1,n1,m2,n2)
%SPHEXP_ACCESS yields array elements of spherical expansion coefficients
%
%   Usage: a = sphexp_access(Anm,m1,n1,m2,n2)
%
%   Input parameters:
%       A           - 1D, 2D, 3D array of expansion coefficients (3rd
%                     dimension for temporal frequency)
%       m1          - order of 1st dimension
%       n1          - degree of 1st dimension (optional, default = |m1|)
%       m2          - order of 2st dimension (optional, default = 0)
%       n2          - degree of 2st dimension (optional, default = |m2|)
%
%   Output parameters:
%       a           - coefficients
%
%   SPHEXP_ACCESS(Anm,m1,n1,m2,n2) computes one/two indices for accessing
%   1D/2D arrays of spherical expansion coefficients
%   A(n1,m1) => A(l1) ; A(n2,m2) => A(l2)
%   CAUTION: THIS FUNCTION DOES NOT USE ANY CHECK OF INPUT ARGUMENTS
%
%   see also: sphexp_access

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

%% ===== Checking of input parameters ===================================
L1 = length(m1);
if nargin < 3
    n1 = abs(m1);
elseif length(n1) < L1
    n1 = repmat(n1,[1 L1]);
elseif length(n1) > L1
    m1 = repmat(m1,[1 length(n2)]);
    L1 = length(m2);
end

if nargin < 4
    m2 = 0;
end

L2 = length(m2);
if nargin < 5
    n2 = abs(m2);
elseif length(n2) < L2
    n2 = repmat(n2,[1 L2]);
elseif length(n2) > L2
    m2 = repmat(m2,[1 length(n2)]);
    L2 = length(m2);
end

%% ===== Computation ====================================================
a = zeros(L1,L2,size(Anm,3));

s1 = abs(m1) <= n1;
s2 = abs(m2) <= n2;

if any(s1) && any(s2)
  [v,w] = sphexp_index(m1(s1),n1(s1),m2(s2),n2(s2));
  a(s1,s2,:) = Anm(v,w,:);
end

end
