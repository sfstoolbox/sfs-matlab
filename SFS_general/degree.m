function phi = degree(phi,conf)
%DEGREE returns the given angle in degree
%
%   Usage: phi = degree(phi,conf)
%          phi = degree(phi)
%
%   Input options:
%       phi     - angle, can be a scalar or matrix (rad)
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output options:
%       phi     - angle (degree)
%
%   DEGREE(phi) returns the given angles phi in degree.
%
%   see also: rad
%

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
phi = phi./pi*180;
