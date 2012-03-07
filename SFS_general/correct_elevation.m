function delta = correct_elevation(delta,conf)
%CORRECT_ELEVATION ensures correct values for elevation angles
%
%   Usage: delta = correct_elevation(delta,conf)
%          delta = correct_elevation(delta)
%
%   Input parameters:
%       delta     - elevation (rad). Can be a single value or a matrix.
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output paramteres:
%       delta     - angle between -pi/2 and +pi/2
%
%   CORRECT_ELEVATION(delta) returns a value for elevation delta between 
%   -pi/2 and pi/2.
%
%   see also: correct_azimuth, read_irs, get_ir
%

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
if nargin==nargmax-1
    conf = SFS_config;
end
if conf.debug
    isargmatrix(delta);
end


%% ===== Computation ====================================================

% Ensure -pi <= delta <= pi
delta = correct_azimuth(delta);

% Ensure -pi/2 <= delta <= pi/2
delta(delta<-pi/2) = -delta(delta<-pi/2) - pi;
delta(delta>pi/2) = -delta(delta>pi/2) + pi;
