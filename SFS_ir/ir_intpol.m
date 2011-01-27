function ir = ir_intpol(ir1,ir2,alpha,conf)
%IR_INTPOL interpolates two given IRs for the given angle
%   Usage: ir = ir_intpol(ir1,ir2,alpha,conf)
%          ir = ir_intpol(ir1,ir2,alpha)
%
%   Input parameters:
%       ir1     - IR with lower abs(angle)
%       ir2     - IR with bigger abs(angle)
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

if nargchk(3,4,nargin)
    error('Wrong number of args. Usage: ir = ir_intpol(ir1,ir2,alpha,conf)');
end

if ~isstruct(ir1)
    error('%s: ir1 has to be a struct!',upper(mfilename));
end
if ~isstruct(ir2)
    error('%s: ir2 has to be a struct!',upper(mfilename));
end
if ~isnumeric(alpha) || ~isscalar(alpha)
    error('%s: alpha has to be a scalar!',upper(mfilename));
end
if nargin<4
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

% Check if the given IR have the same azimuth or the same elevation angle
% in order to get the right interpolation dimension
if ir1.angle(1)==ir2.angle(1) && ir1.angle(2)==ir2.angle(2)
    error('%s: The angles of the two given IRs are the same!',upper(mfilename));
elseif ir1.angle(1)~=ir2.angle(2) && ir1.angle(2)~=ir2.angle(2)
    error(['%s: Both, azimuth and elevation angle are different, ',...
          'no linear interpolation feasible!'],upper(mfilename));
elseif ir1.angle(1)~=ir2.angle(1)  
    
    % Check if the given angle lies between the two IR angles
    % FIXME: this doesn't work easy with the 0..2pi circle!
    %if ~(ir1.angle(1)<alpha && alpha<ir2.angle(1))
    %    error('%s: abs(alpha) has to be > %.3f and < %.3f!',upper(mfilename),...
    %        abs(ir1.angle(1)),abs(ir2.angle(1)));
    %end

    % Interpolation of the azimuth angle
    ir(:,1) = ir1.left + (ir2.left-ir1.left) / ...
                (ir2.angle(1)-ir1.angle(1))*(alpha-ir1.angle(1));
    ir(:,2) = ir1.right + (ir2.right-ir1.right) / ...
                (ir2.angle(1)-ir1.angle(1))*(alpha-ir1.angle(1));
            
elseif ir1.angle(2)~=ir2.angle(2)
    
    % Check if the given angle lies between the two IR angles
    if ~(abs(ir1.angle(2))<abs(alpha)<abs(ir2.angle(2)))
        error('%s: abs(alpha) has to be > %.1f and < %.1f!',upper(mfilename),...
            abs(ir1.angle(2)),abs(ir2.angle(2)));
    end
    
    % Interpolation of the elevation angle
    ir(:,1) = ir1.left + (ir2.left-ir1.left) / ...
                (ir2.angle(2)-ir1.angle(2))*(alpha-ir1.angle(2));
    ir(:,2) = ir1.right + (ir2.right-irs1.right) / ...
                (ir2.angle(2)-ir1.angle(2))*(alpha-ir1.angle(2));
            
end   


%% ===== Plotting ========================================================
% Plot interpolation results (if enabled in SFS_config.m)
if(useplot)
    figure;
    t=1:length(ir1.left);
    plot(t,ir1.left,t,ir(:,1),t,ir2.left);
    legend('IR1','interpolated IR','IR2');
end
