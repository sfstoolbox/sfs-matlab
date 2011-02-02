function ir = ir_intpol(ir1,beta1,ir2,beta2,alpha,conf)
%IR_INTPOL interpolates two given IRs for the given angle
%   Usage: ir = ir_intpol(ir1,beta1,ir2,beta2,alpha,conf)
%          ir = ir_intpol(ir1,beta1,ir2,beta2,alpha)
%
%   Input parameters:
%       ir1     - IR with lower angle
%       beta1   - angle of ir1 (rad)
%       ir2     - IR with bigger angle
%       beta2   - angle of ir2 (rad)
%       alpha   - angle of the desired IR (rad)
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameters:
%       ir      - IR for the given angle alpha (length(IR1),2)
%
%   IR_INTPOL(ir1,ir2,alpha) interpolates the two given IRs ir1 and ir2 for the
%   given angle alpha and returns an interpolated IR.
%
%   see also: get_ir, read_irs
%

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 6;
error(nargchk(nargmin,nargmax,nargin));

if ~isnumeric(ir1) || size(ir1,2)~=2
    error('%s: ir1 has to be a samplesx2 matrix!',upper(mfilename));
end
if ~isnumeric(beta1) || ~isscalar(beta1)
    error('%s: beta1 has to be a scalar!',upper(mfilename));
end
if ~isnumeric(ir2) || size(ir2,2)~=2
    error('%s: ir2 has to be a samplesx2 matrix!',upper(mfilename));
end
if ~isnumeric(beta2) || ~isscalar(beta2)
    error('%s: beta2 has to be a scalar!',upper(mfilename));
end
if ~isnumeric(alpha) || ~isscalar(alpha)
    error('%s: alpha has to be a scalar!',upper(mfilename));
end
if nargin<nargmax
    useconfig = true;
elseif ~isstruct(conf)
    error('%s: conf has to be a struct.',upper(mfilename));
else
    useconfig = false;
end


%% ===== Configuration ==================================================

% Load default configuration values
if(useconfig)
    conf = SFS_config;
end

useplot = conf.useplot;


%% ===== Computation ====================================================

% === HRIR interpolation ===

% IR interpolation for left and right signal after the following
% formula (_s means smaller and _l larger):
% int_hrir = hrir_s + (hrir_l-hrir_s) / (phi_l-phi_s) * ...
%   (phi-phi_s)

% Check if the given IR have the same angle
% in order to get the right interpolation dimension
if beta1==beta2
    error('%s: The angles of the two given IRs are the same!',upper(mfilename));
end

% Correct the given angles
beta1 = correct_azimuth(beta1);
beta2 = correct_azimuth(beta2);
alpha = correct_azimuth(alpha);

if alpha==beta1 || alpha==beta2
    error('%s: no interpolation needed for the given alpha value.',...
        upper(mfilename));
end
if (alpha>beta1 && alpha>beta2) || (alpha<beta1 && alpha<beta2)
    error('%s: alpha has to be between beta1 and beta2.',upper(mfilename));
end
if length(ir1)~=length(ir2)
    error('%s: the given IRs have not the same length.',upper(mfilename));
end

% Look for the angle where alpha is nearest to and use that IR as the starting
% point for the interpolation
if abs(alpha-beta1)<=abs(alpha-beta2)

    ir = ir1 + (ir2-ir1) ./ (beta2-beta1)*(alpha-beta1);

else

    ir = ir2 + (ir1-ir2) ./ (beta1-beta2)*(alpha-beta2);

end


%% ===== Plotting ========================================================
% Plot interpolation results (if enabled in SFS_config.m)
if(useplot)
    figure;
    t=1:length(ir);
    plot(t,ir1(:,1),t,ir(:,1),t,ir2(:,1));
    legend('IR1','interpolated IR','IR2');
end
