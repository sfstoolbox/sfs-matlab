function out = sphbessely(nu,z)
% SPHBESSELY spherical bessel function of second kind, of order nu, and argument z
%   Usage: out = sphbesselj(nu,z)
%
%   Input parameters:
%       nu  - order of bessel function
%       z   - argument of bessel function
%
%   Output parameters:
%       out - value of bessel function at point z
%
%   SPHBESSELY(nu,z) spherical bessel function of order nu, frist type, and 
%   argument z
%
%   see also: sphbesselh, sphbesselj
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
out = sqrt(pi./(2.*z)) .* bessely(nu+0.5, z);
