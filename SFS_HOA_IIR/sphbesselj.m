function [ out ] = sphbesselj(nu, z)
%function [ out ] = sphbesselj(nu, z)
% SPHBESSELJ spherical bessel function of first type of order nu and 
% argument z
% Author: J. Ahrens, 09.04.2008

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
