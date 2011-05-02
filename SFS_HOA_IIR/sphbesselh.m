function [ out ] = sphbesselh(nu, k, z)
% SPHBESSELH spherical hankel function of order nu, kind k, and argument z
%   Usage: [ out ] = sphbesselh(nu, k, z)
%
%   Input parameters:
%       nu  - order of bessel function
%       k   - kind of bessel function (1 ^= first kind, 2 ^= second kind)
%       z   - argument of bessel function
%
%   Output parameters:
%       out - value of bessel function at point z
%
%   SPHBESSELH(nu,k,z) spherical hankel function of order nu, kind k, and 
%   argument z
%
%   see also: sphbesselj, sphbessely
%

% Author: J. Ahrens, 09.04.2008

% ------ Checking of input parameters ---------
  
if nargchk(3,3,nargin)
    error(['Wrong number of args. ',...
           'Usage: [ out ] = sphbesselh(nu, k, z)']);
end

if ~isnumeric(nu) || ~isscalar(nu) || nu<0
    error('%s: nu must be a positive scalar.',upper(mfilename));
end;

if (k==1) 
    sign = 1;
elseif (k==2)
    sign = -1;
else
    error(['%s: Invalid kind of Hankel function is asked ',...
           '(k has to be 1 or 2).'],upper(mfilename));
end

if ~isnumeric(z)
    error('%s: z has to be numeric.');
end


% ------ Computation -------------------------

out = sphbesselj(nu, z) + 1j .* sign .* sphbessely(nu, z);
