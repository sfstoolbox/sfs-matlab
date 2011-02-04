function phi = degree(phi)
%DEGREE returns the given angle in degree
%   Usage: phi = degree(phi)
%
%   Input options:
%       phi - angle (rad)
%
%   Output options:
%       phi - angle (degree)
%
%   DEGREE(phi) returns the given angle phi in degree.
%
%   see also: rad

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
if ~isnumeric(phi) || ~isscalar(phi)
    error('%s: phi must be a scalar.',upper(mfilename));
end


%% ===== Computation =====================================================
phi = phi/pi*180;
