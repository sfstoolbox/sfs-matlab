function delta = correct_elevation(delta)
%CORRECT_ELEVATION ensures correct values for elevation angles 
%   Usage: delta = correct_elevation(delta)
%
%   Input parameters:
%       delta     - elevation (rad). Can be a single value or a vector.
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


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
isargvector({delta},{'delta'});


%% ===== Computation ====================================================

% Ensure -pi <= delta <= pi
delta = correct_azimuth(delta);

% Ensure -pi/2 <= delta <= pi/2
delta(delta<-pi/2) = -delta(delta<-pi/2) - pi;
delta(delta>pi/2) = -delta(delta>pi/2) + pi;
