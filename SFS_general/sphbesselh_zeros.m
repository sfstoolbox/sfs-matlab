function [z,p] = sphbesselh_zeros(order)
%SPHBESSELH_ZEROS finds zeros/roots of spherical hankel function
%
%   Usage: [z,p] = sphbesselh_zeros(order)
%
%   Input parameters:
%       order       - order of hankel function
%
%   Output parameters:
%       z       - zeros/roots fo the Bessel function
%       p       - roots of the Bessel function
%
%   SPHBESSELH_ZEROS(order) finds zeros and roots for a spherical hankel function
%   of the specified order. Due to numerical problems, the order is limited up
%   to 85.
%
%   See also: sphbesselh, driving_function_imp_nfchoa

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


%% ===== Checking input parameters =======================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
isargpositivescalar(order);


%% ===== Configuration ===================================================
% Method for calculating zeros
% See https://github.com/sfstoolbox/sfs/issues/57 for a discussion
method = 1;


%% ===== Main ============================================================
if order<86
    if method==1
        % Method 1 after FIXME
        % Formula for nominator (source?)
        B = zeros(1,order+2);
        for n=0:order
            B(n+1) = factorial(2*order-n)/(factorial(order-n)*factorial(n)*2^(order-n));
        end
        B = B(end:-1:1);
        % Zeros
        z = roots(B);
        % Poles (are always zero)
        p = zeros(order,1);
    elseif method==2
        % Method 2 after Pomberger 2008, p. 43
        B = cell(order+1,1);
        z = cell(order+1,1);
        p = cell(order+1,1);
        for n=0:order
          % Recursion formula for nominator
          B{n+1} = zeros(1,n+1);
          for k=0:n-1
              B{n+1}(k+1) = ((2*n-k-1)*(2*n-k)) / (2*(n-k)) * B{n}(k+1);
          end
          B{n+1}(n+1) = 1;
        end
        for n=0:order
          % Flip nominator polynoms
          B{n+1} = B{n+1}(end:-1:1);
          % Zeros
          z{n+1} = roots(B{n+1});
          % Poles (are always zero)
          p{n+1} = zeros(order,1);
        end
        z = z{order+1};
        p = p{order+1};
    else
        error('%s: method has to be 1 or 2.',upper(mfilename));
    end
else
    error(['%s: for orders higher than 85 no stable numerical ', ...
           'method is available at the moment to caclulate the zeros.'], ...
          upper(mfilename));
end
return


%% ===== Computation with Multiprecission Toolbox ========================
% For the Multiprecission Toolbox, see: http://www.advanpix.com
% Unfortunately it turned out, that the obtained zeros with this method have
% some systematic errors, see
% https://github.com/sfstoolbox/sfs/issues/57#issuecomment-183791477
% The following code was used to calculate the zeros with the Multiprecission
% Toolbox. The results are stored at
% https://github.com/sfstoolbox/data/tree/master/sphbesselh_zeros
B = mp(zeros(1,order+2));
A = B;
for n=mp(0:order)
    B(n+1) = factorial(2*order-n)/(factorial(order-n)*factorial(n)*2^(order-n));
end
B = B(end:-1:1);
z = roots(B);
A(2) = mp(1);
p = roots(A);
