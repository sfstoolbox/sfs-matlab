function [z,p] = sphbesselh_zeros(order)
%SPHBESSELH_ZEROS finds zeros/roots of spherical hankel function
%
%   Usage: [z,p] = sphbesselh_zeros(order)
%
%   Input parameters:
%       order   - order of hankel function
%
%   Output parameters:
%       z       - zeros/roots fo the Bessel function
%       p       - roots of the Bessel function
%
%   SPHBESSELH_ZEROS(order) finds zeros and poles for a spherical hankel 
%   function of the specified order. It is based on the investigations in 
%   Hahn and Spors (2017) and the Python implementation in from scipy in
%   signal.filter_design._bessel_zeros.
%
%   See also: sphbesselh, driving_function_imp_nfchoa
%
%   References:
%
%       N. Hahn, S. Spors (2017) - "Further Investigations on the Design of 
%       Radial Filters for the Driving Functions of Near-Field Compensated 
%       Higher-Order Ambisonics," in 142nd Convention of the Audio Engineering
%       Society, Paper 9732, http://www.aes.org/e-lib/browse.cfm?elib=18609
%
%       R. Campos, M. Calderon (2011) - "Approximate closed-form formulas for
%       the zeros of the Bessel Polynomials," https://arxiv.org/abs/1105.0957
%
%       This implementation is based on scipy: http://bit.ly/2tPfePn

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
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
isargpositivescalar(order);


%% ===== Constants =======================================================
TOL = 1E-15;
MAXITER = 50;


%% ===== Main ============================================================
p = zeros(order,1);  % poles are always zero

if order < 2
    z = -ones(order,1);
    return
end

% Approximate roots of nth-order ordinary Bessel polynomial as starting points
r = campos_zeros(order);

% Zeros of nth-order ordinary Bessel polynomial y_n are the same as for
% the exponentially scaled modified Bessel function of second kind K_v:
% \sqrt{2pi/x} * exp(1/x) * K_{n+0.5}(1/x) = y_n(x),
% see for example http://bit.ly/2unx63m.
% Hence, we define the target function and its first derivative as
f  = @(x) besselk(order+0.5, 1./x, 1);
fp = @(x) besselk(order-0.5, 1./x, 1)./2./x.^2  - ...
          besselk(order+0.5, 1./x, 1)./x.^2  + ...
          besselk(order+1.5, 1./x, 1)./2./x.^2;

% Simulataneous root finding using the method of Aberth-Ehrlich method
r = aberth(f,fp,r,TOL,MAXITER);

% Refine root position using Newton-Raphson method
for idx=1:order
    r(idx) = newton(f,fp,r(idx),TOL,MAXITER);
end

% Average complex conjugates to make them exactly symmetrical
r = 0.5*r + 0.5*conj(r(end:-1:1));

% Roots should sum to -1
if abs(sum(r) + 1) > TOL
    error(['%s: Generated roots of the ordinary Bessel polynomial are' ...
      ' inaccurate, order=%d'],upper(mfilename),order);
end

% Zeros are the inverted roots of nth-order ordinary Bessel polynomial
z = 1./r;

end


%% ===== Auxiliary Functions =============================================
function z0 = campos_zeros(n)
    % Approximate roots of ordinary Bessel polynomial of nth order, see
    % Campos and Calderon (2011) and http://bit.ly/2txLJyM

    if n == 1
        z0 = -1;
        return
    end

    % See Campos and Calderon (2011), Eq. (7)
    r = polyval( [1 2 0 0],n);          % n^3 + n^2
    a1 = polyval( [-6 -6],n) / r;       % ( -6n - 6 ) / r
    a2 = 6 / r;

    s = polyval( [1 -3 0 2 0 0],n);     % n^5 - 3n^4 + 2n^2
    b3 = polyval( [-8 16],n) / s;       % ( -8n + 16 ) / s
    b2 = polyval( [12 -12 -24],n) / s;  % ( 12n^2 - 12n - 24 ) / s
    b1 = polyval( [-2 -12 24 8],n) / s; % ( -2n^3 - 12n^2 + 24n + 8 ) / s
    b0 = polyval( [-1 5 0 -6 0],n) / s; % ( -n^4 + 5n^3 - 6n ) / s

    % See Campos and Calderon (2011), Eq. (4)
    k = (1:n).';
    x0 = polyval( [a2 a1 0],k);         % real part
    y0 = polyval( [b3 b2 b1 b0],k);     % imaginary part

    % See Campos and Calderon (2011), Eq. (8)
    z0 = x0 + 1j*y0;

end

function x = aberth(f,fp,x0,TOL,MAXITER)
    % Ehrlich-Aberth method to simulatenous approximation of roots,
    % see http://bit.ly/2sOfeiT
    N = length(x0);
    beta = zeros(size(x0));
    x = x0;
    for iter=1:MAXITER
        alpha = -f(x) ./ fp(x);
        
        for k=1:N
            beta(k) = sum(1./(x(k) - x([1:k-1,k+1:N])));
        end
        
        x = x + alpha./(1 + alpha .* beta);
        
        if all(abs(alpha) <= TOL)
            return
        end
    end
    warning('%s: maximum iterations in aberth() reached; order = %d', ... 
        upper(mfilename), N);
end

function x = newton(f,fp,x0,TOL,MAXITER)
    % Newton-Rapheson method for approximation of a single root
    for iter=1:MAXITER
       x = x0 - f(x0) ./ fp(x0);   
       if abs(x - x0) < TOL
           return
       end
       x0 = x;   
    end
    warning('%s: maximum iterations in newton() reached',  upper(mfilename));
end
