function [b, a] = thiran_filter(Norder, fdt)
%THIRAN_FILTER computes Thiran's IIR allpass for Maximally Flat Group Delay
%
%   Usage: [b, a] = lagrange_filter(Norder, fdt)
%
%   Input parameter:
%     Norder - order of filter
%     fdt    - vector of fractional delays -0.5 <= fdt < 0.5
%
%   Output parameter:
%     b   - numerator polynomial of H(z) / [Norder+1 x Nfdt]
%     a   - denominator polynomial of H(z) / [Norder+1 x Nfdt]
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

% shift fractional delay in order to optimize performance
fdt = fdt(:);  % ensure column vector
Nfdt = numel(fdt);

% denomimator polynomial of H(z)
a = [ones(1,Nfdt); zeros(Norder, Nfdt)];
for kdx=1:Norder
  a(kdx+1,:) = (-1).^kdx * ...
    factorial(Norder)/(factorial(kdx)*factorial(Norder-kdx)) * ...
    prod( bsxfun(@plus, fdt, 0:Norder)./bsxfun(@plus, fdt, kdx:kdx+Norder), 2 );        
end
% numerator polynomial of H(z)
b = a(end:-1:1,:);

end
