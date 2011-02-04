function phi = degree(phi)
%DEGREE returns the given angle in degree
%   Usage: phi = degree(phi)
%
%   Input options:
%       phi - angle, can be a scalar or matrix (rad)
%
%   Output options:
%       phi - angle (degree)
%
%   DEGREE(phi) returns the given angles phi in degree.
%
%   see also: rad

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
if ~isnumeric(phi) || ~ismatrix(phi)
    error('%s: phi must be a scalar.',upper(mfilename));
end


%% ===== Computation =====================================================
phi = phi./pi*180;
