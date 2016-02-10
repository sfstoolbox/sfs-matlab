function [l1, l2] = sphexp_index(m1, n1, m2, n2)
%SPHEXP_INDEX calculates index(indices) for the arrays of spherical expansion 
%coefficients
%
%   Usage: [l1, l2] = sphexp_index(m1, n1, m2, n2)
%
%   Input parameters:
%       m1          - order of 1st dimension
%       n1          - degree of 1st dimension (optional, default = |m1|)
%       m2          - order of 2st dimension (optional, default = 0)
%       n2          - degree of 2st dimension (optional, default = |m2|)
%
%   Output parameters:
%       l1          - index for 1st dimension
%       l2          - index for 2st dimension
%
%   SPHEXP_INDEX(m1, n1, m2, n2) computes one/two indices for accessing
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

%% ===== Checking of input  parameters ==================================
if nargin < 2
  n1 = abs(m1);
end
if nargin < 3
  m2 = 0;
end
if nargin < 4 
  n2 = abs(m2);
end

%% ===== Computation ====================================================
l1 = (n1 + 1).^2 - (n1 - m1);
l2 = (n2 + 1).^2 - (n2 - m2);
