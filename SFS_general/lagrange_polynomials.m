function [pi, L] = lagrange_polynomials(xi, yi)
%LAGRANGE_POLYNOMIALS compute lagrange polynomials based on sampling positions
%
%   Usage: [li, [L]] = lagrange_polynomials(xi, [yi])
%
%   Input parameter:
%     xi - row vector sampling positions x_i with i = 0, ..., N
%     yi - optional row vector of values of sampling positions y_i = f(x_i)
%
%   Output parameter:
%     pi  - matrix of lagrange polynomials p_i(x) with i = 0, ..., N
%           p_i(x) = (x - x_0)/(x_i - x_0) * ... * (x - x_i-1) /(x_i - x_i-1) *
%                    (x - x_i+1) /(x_i - x_i+1) ... (x - x_N)/(x_i - x_N)
%                  = c_{i,N} x^N + c_{i, N-1} x^(N-1) + ... + c_{i, 0} x^(0)
%           each row of pi stores the coefficients c_{i,k} with k = N, ..., 0 
%     L   - optional interpolation polynomial L(x)
%           L(x) = y_0 * p_0(x) + y_1 * p_1(x) + ... + y_N * p_N(x)
%
%   LAGRANGE_POLYNOMIALS computes N interpolation polynomials, each of N-th
%   order, based on N sampling positions.
%
%   See also: delayline_read, delayline_write

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

%% ===== Computation =====================================================

% number elements in x_i determines degree of the lagrange polynoms
N = length(xi);  

% each row contains [1 x_i] which is equivalent to m_i(x) = (x - x_i)
mi = [ones(N,1), -xi'];

% l(x) = m_0(x) * m_1(x) ... * m_N(x) = (x - x_0) * (x - x_1) * ... * (x - x_N)
l = 1;
for idx=1:N
  % convolution of coefficients means multiplication of polynoms
  l = conv(l,mi(idx,:));  
end

pi = zeros(N,N);
for idx=1:N
  % nom_i(x) = l(x) / m_i(x) 
  %          = (x - x_0) * ... * (x - x_{i-1}) * (x - x_{i+1}) * ... * (x - x_N)  
  nominator = deconv(l,mi(idx,:));
  % denom_i = nom_i(x_i) evaluated with "polyval"
  denominator = polyval(nominator,xi(idx));  
  % l_i(x) = nom_i(x) / denom_i;
  pi(idx,:) = nominator./denominator;
end

% optional computation of interpolation polynom L(x)
% L(x) = y_0 * l_0(x) + y_1 * l_1(x) + ... + y_N * l_N(x)
if nargin == 2 && nargout == 2
  L = yi * pi;  % row vector * matrix = row vector
end