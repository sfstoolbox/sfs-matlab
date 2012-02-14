function out = sphbesselj(nu,z)
% SPHBESSELJ spherical bessel function of first kind of order nu, and argument z
%   Usage: out = sphbesselj(nu,z)
%
%   Input parameters:
%       nu  - order of bessel function
%       z   - argument of bessel function
%
%   Output parameters:
%       out - value of bessel function at point z
%
%   SPHBESSELJ(nu,z) spherical bessel function of order nu, frist type, and 
%   argument z
%
%   see also: sphbesselh, sphbessely
%
% AUTHOR: Jens Ahrens
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking input parameters =======================================
nargmin = 2;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargpositivescalar(nu)
isargnumeric(z)


%% ===== Computation =====================================================
out = zeros(size(z));

% avoid division by "0"
if (nu==0)
    out(find(z==0)) = 1;
elseif (nu~=0)
    out(find(z==0)) = 0;
end

% finally evaluate for z~=0
out(find(z~=0)) = sqrt(pi./(2.*z(find(z~=0)))) .* ...
                                    besselj(nu+0.5, z(find(z~=0)));
