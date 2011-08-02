function out = sphbesselh(nu,k,z)
% SPHBESSELH spherical hankel function of order nu, kind k, and argument z
%   Usage: out = sphbesselh(nu,k,z)
%
%   Input parameters:
%       nu  - order of hankel function
%       k   - kind of hankel function (1 ^= first kind, 2 ^= second kind)
%       z   - argument of hankel function
%
%   Output parameters:
%       out - value of hankel function at point z
%
%   SPHBESSELH(nu,k,z) spherical hankel function of order nu, kind k, and 
%   argument z
%
%   see also: sphbesselj, sphbessely

% AUTHOR: Jens Ahrens


%% ===== Checking input parameters =======================================
nargmin = 3;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));
isargpositivescalar(nu)
if (k==1)
    sign = 1;
elseif (k==2)
    sign = -1;
else
    error(['%s: Invalid kind of Hankel function is asked ',...
           '(k has to be 1 or 2).'],upper(mfilename));
end
isargnumeric(z)


%% ===== Computation =====================================================
out = sphbesselj(nu, z) + 1j .* sign .* sphbessely(nu, z);
