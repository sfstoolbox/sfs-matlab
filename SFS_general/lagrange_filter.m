function [b, li] = lagrange_filter(Norder, fdt)
%LAGRANGE_FILTER computes Lagrange interpolation filter for fractional delays
%
%   Usage: [b, [li]] = lagrange_filter(Norder, [fdt])
%
%   Input parameter:
%     Norder - order of Lagrange polynomials
%     fdt    - optional vector of fractional delays
%              0 <= fdt < 1 if Norder is odd,
%              -0.5 <= fdt < 0.5 if Norder is even,
%
%   Output parameter:
%     b   - filter coefficients / [Norder+1 x Nfdt]
%     li  - optional matrix of lagrange polynomials l_i(x) with i = 0, ..., N
%           l_i(x) = (x - x_0)/(x_i - x_0) * ... * (x - x_i-1) /(x_i - x_i-1) *
%                    (x - x_i+1) /(x_i - x_i+1) ... (x - x_N)/(x_i - x_N)
%                  = c_{i,N} x^N + c_{i, N-1} x^(N-1) + ... + c_{i, 0} x^(0)
%           where x_k = k with k = N, ..., 0. The i-th row of li stores the 
%           coefficients c_{i,k}.
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

if nargin == 2
  D = fdt(:).' + floor(Norder/2);
  
  % aux = D, (D-1), (D-2),...,(D-N+1),(D-N)
  aux = bsxfun(@minus, D, (0:Norder).');  
  
  % denom = n*(n-1)*...*(n-N+1)*(n-N) = n!*(N-n)!*(-1)^(N-n)
  denom = factorial(0:Norder);
  denom = denom.*denom(end:-1:1).*(-1).^(Norder:-1:0);
  
  b = zeros(Norder+1, length(D));
  for ndx=1:Norder+1
    b(ndx,:) = prod(aux([1:ndx-1,ndx+1:end],:), 1)./denom(ndx);
  end
else
  b = [];
end

if nargout == 2  
  % each row contains [1 x_i] which is equivalent to m_i(x) = (x - x_i)
  mi = [ones(Norder,1), -(0:Norder).'];
  
  % l(x) = m_0(x) * m_1(x) .. * m_N(x) = (x - x_0) * (x - x_1) * .. * (x - x_N)
  l = 1;
  for idx=1:Norder
    % convolution of coefficients means multiplication of polynoms
    l = conv(l,mi(idx,:));
  end
  
  li = zeros(Norder,Norder);
  for idx=1:N
    % nom_i(x) = l(x) / m_i(x)
    %          = (x - x_0) * .. * (x - x_{i-1}) * (x - x_{i+1}) * .. * (x - x_N)
    nominator = deconv(l,mi(idx,:));
    % denom_i = nom_i(x_i) evaluated with "polyval"
    denominator = polyval(nominator,xi(idx));
    % l_i(x) = nom_i(x) / denom_i;
    li(idx,:) = nominator./denominator;
  end   
end
