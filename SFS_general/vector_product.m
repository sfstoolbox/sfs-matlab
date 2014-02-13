function y = vector_product(x1,x2,dim)
%VECTOR_PRODUCT calculates the scalar product for two vectors within matrices
%
%   Usage: y = vector_product(x1,x2,dim)
%
%   Input parameters:
%       x1  - first matrix [n x m]
%       x2  - second matrix [n x m]
%       dim - dimension along the scalar product should be calculated
%
%   Output parameter:
%       y   - scalar product between the vectors [1 x m] or [n x 1]
%
%   VECTOR_PRODUCT(x1,x2,dim) calculates the scalar product between the vectors
%   given within the matrices along the dimension dim.
%
%   see also: vector_norm

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
nargmin = 3;
nargmax = 3;
narginchk(nargmin,nargmax)
isargmatrix(x1,x2)
isargpositivescalar(dim)
if size(x1)~=size(x2)
    error('%s: the given matrices x1, x2 had to be of the same size.', ...
        upper(mfilename));
end


%% ===== Computation =====================================================
y = sum(x1.*x2,dim);
