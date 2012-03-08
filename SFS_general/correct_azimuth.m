function phi = correct_azimuth(phi,conf)
%CORRECT_AZIMUTH ensures correct values for azimuth angles
%
%   Usage: phi = correct_azimuth(phi,conf)
%          phi = correct_azimuth(phi)
%
%   Input parameters:
%       phi     - azimuth (rad). Can be a single value or a matrix.
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output paramteres:
%       phi     - angle between -pi and +pi-eps
%
%   CORRECT_AZIMUTH(phi) returns a value for azimuth phi between 
%   -pi and +pi-eps.
%
%   see also: correct_elevation, read_irs, get_ir
%

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
if nargin==(nargmax-1)
    conf = SFS_config;
end
if conf.debug
    isargmatrix(phi);
end


%% ===== Computation ====================================================

% Ensure -2pi <= phi <= 2pi
phi = rem(phi,2*pi);

% FIXME: I think for the OpenDAFF format the opposite is needed for pi and -pi!

% Ensure -pi <= phi < pi
phi(phi<-pi) = phi(phi<-pi) + 2*pi;
phi(phi>=pi) = phi(phi>=pi) - 2*pi;
