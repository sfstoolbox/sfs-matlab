function y = vector_norm(x,dim)
%VECTOR_NORM calculates the norm for the rows/columns of a matrix
%
%   Usage: y = vector_norm(x,dim)
%
%   Input parameters:
%       x   - matrix [n x m]
%       dim - dimension along the norm should be calculated
%
%   Output parameter:
%       y   - norm of the vectors [1 x m] or [n x 1]
%
%   VECTOR_NORM(x,dim) calculates the p-norm (with p=2) for the vectors
%   given within the matrix along the dimension dim.
%
%   see also: norm, vector_product

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
% NOTE: this is disabled due to performance issues in HRTF extrapolation, where
% this function is called multiple times.
%nargmin = 2;
%nargmax = 2;
%narginchk(nargmin,nargmax)
%isargmatrix(x)
%isargpositivescalar(dim)


%% ===== Computation =====================================================
y = sum(abs(x).^2,dim).^(1/2);
