function [phi,delta] = correct_angle(phi)
%CORRECT_ANGLE Ensures correct values for azimuth angles 
%   Usage: phi = correct_angle(phi)
%
%   Input parameters:
%       phi     - azimuth (rad)
%
%   Output paramteres:
%       phi     - angle between -pi and pi-eps
%
%   CORRECT_ANGLE(phi) returns a value for azimuth phi between -pi and pi-eps.
%
%   see also: read_irs, get_ir
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================

if nargchk(1,1,nargin)
    error('Wrong number of args. Usage: phi = correct_azimuth(phi)');
end
if ~isnumeric(phi) || ~isvector(phi)
    error('%s: phi has to be a vector!',upper(mfilename));
end


%% ===== Computation ====================================================

% Ensure -2pi <= phi <= 2pi
phi = rem(phi,2*pi);

% NOTE: the choice of -pi <= phi < pi has been taken due to the same in the
% HRIR data set of Alexander Lindau and because we start at the lowest
% value -pi.
% FIXME: I think for the OpenDAFF format the opposite is needed for pi and -pi!

% Ensure -pi <= phi < pi
phi(phi<-pi) = phi(phi<-pi) + 2*pi;
phi(phi>=pi) = phi(phi>=pi) - 2*pi;
