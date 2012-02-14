function ir = intpol_ir(ir1,beta1,ir2,beta2,alpha)
%INTPOL_IR interpolates two given IRs for the given angle
%   Usage: ir = intpol_ir(ir1,beta1,ir2,beta2,alpha)
%
%   Input parameters:
%       ir1     - IR with lower angle
%       beta1   - angle of ir1 (rad)
%       ir2     - IR with bigger angle
%       beta2   - angle of ir2 (rad)
%       alpha   - angle of the desired IR (rad)
%
%   Output parameters:
%       ir      - IR for the given angle alpha (length(IR1),2)
%
%   INTPOL_IR(ir1,beta1,ir2,beta2,alpha) interpolates the two given IRs ir1 and 
%   ir2 with their corresponding angles beta1 and beta2 for the given angle 
%   alpha and returns an interpolated IR.
%
%   see also: get_ir, shorten_ir, read_irs
%

% AUTHOR: Sascha Spors, Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 5;
error(nargchk(nargmin,nargmax,nargin));

isargscalar(beta1,beta2,alpha);

if ~isnumeric(ir1) || size(ir1,2)~=2
    error('%s: ir1 has to be a samples x 2 matrix!',upper(mfilename));
end
if ~isnumeric(ir2) || size(ir2,2)~=2
    error('%s: ir2 has to be a samples x 2 matrix!',upper(mfilename));
end


%% ===== Computation ====================================================

% Check if the given IR have the same angle
% in order to get the right interpolation dimension
if beta1==beta2
    error('%s: The angles of the two given IRs are the same!',upper(mfilename));
end

% Correct the given angles
beta1 = correct_azimuth(beta1);
beta2 = correct_azimuth(beta2);
alpha = correct_azimuth(alpha);

%if alpha==beta1 || alpha==beta2
%    error('%s: no interpolation needed for the given alpha value.',...
%        upper(mfilename));
%end
if length(ir1)~=length(ir2)
    error('%s: the given IRs have not the same length.',upper(mfilename));
end

% Linear interpolate the two given IRs
ir = ir1 + (ir2-ir1) ./ (beta2-beta1)*(alpha-beta1);
