function [b,a] = thiran_filter(Norder,fdt)
%THIRAN_FILTER computes Thiran's IIR allpass for Maximally Flat Group Delay
%
%   Usage: [b,a] = thiran_filter(order,fdt)
%
%   Input parameter:
%     order  - order of filter
%     fdt    - vector of fractional delays -0.5 <= fdt < 0.5
%
%   Output parameter:
%     b   - numerator polynomial of H(z) / [Norder+1 x Nfdt]
%     a   - denominator polynomial of H(z) / [Norder+1 x Nfdt]
%
%   See also: delayline, lagrange_filter

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Computation =====================================================

% Shift fractional delay in order to optimize performance
fdt = fdt(:);  % ensure column vector
Nfdt = numel(fdt);

% Denomimator polynomial of H(z)
a = [ones(1,Nfdt); zeros(Norder,Nfdt)];
for kdx=1:Norder
    a(kdx+1,:) = (-1).^kdx * ...
        factorial(Norder)/(factorial(kdx)*factorial(Norder-kdx)) * ...
        prod( bsxfun(@plus,fdt,0:Norder)./bsxfun(@plus,fdt,kdx:kdx+Norder), 2 );
end
% Numerator polynomial of H(z)
b = a(end:-1:1,:);
