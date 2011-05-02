function [ out ] = sphbessely(nu, z)
%function [ out ] = sphbessely(nu, z)
%
%SPHBESSELY spherical bessel function of second kind of order nu and 
% argument z
% Author: J. Ahrens, 09.04.2008

out = sqrt(pi./(2.*z)) .* bessely(nu+0.5, z);
