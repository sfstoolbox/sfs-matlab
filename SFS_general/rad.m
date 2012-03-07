function phi = rad(phi,conf)
%RAD returns the given angle in radians
%
%   Usage: phi = rad(phi,conf)
%          phi = rad(phi)
%
%   Input options:
%       phi     - angle, can be a scalar or matrix (degree)
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output options:
%       phi     - angle (rad)
%
%   RAD(phi) returns the given angles phi in radians.
%
%   see also: degree

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
if nargin==nargmax-1
    conf = SFS_config;
end
if conf.debug
    isargmatrix(phi);
end


%% ===== Computation =====================================================
phi = phi./180*pi;
