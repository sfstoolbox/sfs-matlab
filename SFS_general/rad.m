function phi = rad(phi)
%RAD returns the given angle in radians
%   Usage: phi = rad(phi)
%
%   Input options:
%       phi - angle (degree)
%
%   Output options:
%       phi - angle (rad)
%
%   RAD(phi) returns the given angle phi in radians.
%
%   see also: degree

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
if ~isnumeric(phi) || ~isscalar(phi)
    error('%s: phi must be a scalar.',upper(mfilename));
end


%% ===== Computation =====================================================
phi = phi/180*pi;
