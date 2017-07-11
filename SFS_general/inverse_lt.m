function A = inverse_lt(Am,x)
%INVERSE_LT computes the inverse Legendre transform (ILT)
%
%   Usage: [A,mu] = inverse_lt(Am,Nx)
%
%   Input parameters:
%       Am      - Legendre expansion coefficients [N x M]
%       x       - values for which the ILT is computed [-1:1], [1 x Nx]
%
%   Output parameters:
%       A       - ILT corresponding to mu [N x Nx]
%
%   See also: inverse_cht

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Developers                             *
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


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargmatrix(Am);
isargvector(x);

%% ===== Computation ==================================================
M = (size(Am,2)-1)/2;
N = size(Am,1);
Nx = length(x);

x = x(:).';  % ensure row vector

% Implementation of
%         ___
%         \       2m+1
% A(x) =  /__    ------ A  P (x)
%        m=0..M    2     m  m
%
% with P_m(x) being the mth-order Legendre polynomial

A = zeros(N, Nx);
for m=0:M
    A = A + bsxfun(@times, Am(:,m+1), (2.*m+1)./2.*asslegendre(m,0,x));
end
