function [ Lnm ] = asslegendre(n,m,arg)
%ASSLEGENDRE calculates associates Legendre function of degree n and order m
%
%   Usage: Lnm = asslegendre(n,m,arg)
%
%   Input parameters:
%       n     - spherical harmonic degree
%       m     - spherical harmonic order
%       arg   - values to calculate the lengendre polynomes
%
%   Output parameters:
%       Lnm   - associates Legendre function
%
%   ASSLEGENDRE(n,m,arg) calculates the associates Legendre function of degree
%   n and order m for the values arg.
%
%   See also: sphharmonics, legendre

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


%% ===== Checking input parameters =======================================
nargmin = 3;
nargmax = 3;
narginchk(nargmin,nargmax);
isargpositivescalar(n)
isargscalar(m)
isargnumeric(arg)
if n<abs(m)
    warning( 'Absolute value of order m must be less than or equal to the degree n.' ); 
    Lnm = zeros(size(alpha));
    return;
end


%% ==== Computation ======================================================
Lnm = legendre(n,arg);
if n~=0
    Lnm = squeeze(Lnm(abs(m)+1,:,:));  
end
if m<0
    Lnm = -1.^abs(m) .* factorial(n-abs(m)) ./ factorial(n+abs(m)) .* Lnm;
end
Lnm = reshape(Lnm,size(arg));
