function [weight,delay] = driving_function_imp_wfs_25d(x0,xs,src,conf)
%DRIVING_FUNCTION_IMP_WFS_25D calculates the WFS 2.5D weighting and delaying
%
%   Usage: [weight,delay] = driving_function_imp_wfs_25d(x0,xs,src,conf);
%          [weight,delay] = driving_function_imp_wfs_25d(x0,xs,src);
%
%   Input parameters:
%       x0      - position  and direction of secondary source (m)
%       xs      - position of virtual source or diirection of plane wave (m)
%       src     - source type of the virtual source
%                     'pw' - plane wave (xs, ys are the direction of the
%                            plane wave in this case)
%                     'ps' - point source
%                     'fs' - focused source
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       weight  - weight (amplitude) of the driving function
%       delay   - delay of the driving function (s)
%
%   DRIVING_FUNCTION_IMP_WFS_25D(x0,xs,src,conf) returns the
%   weighting and delay parameters of the WFS 2.5D driving function for the given
%   source type and position and loudspeaker positions.
%
%   see also: wave_field_imp_wfs_25d, driving_function_mono_wfs_25d

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 4;
error(nargchk(nargmin,nargmax,nargin));
isargsecondarysource(x0)
isargposition(xs);
xs = position_vector(xs);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
% Sampling rate
fs = conf.fs;
% Speed of sound
c = conf.c;
% Array type
array = conf.array;
xref = position_vector(conf.xref);


%% ===== Computation =====================================================
% Check also the activity of the used loudspeaker.
ls_activity = secondary_source_selection(x0,xs,src);
if ls_activity>0

    % Direction and position of secondary sources
    nx0 = secondary_source_direction(x0);
    x0 = x0(1:3);

    % Constant amplitude factor
    g0 = sqrt(2*pi*norm(xref-x0));
    
    if strcmp('pw',src)
        % === Plane wave ===
        % Direction of plane wave
        nxs = xs / norm(xs); 
        % Delay and amplitude weight
        % NOTE: <n_pw,n(x0)> is the same as the cosinus between their angle
        delay = 1/c * nxs*x0';
        weight = 2*g0 * nxs*nx0';

    elseif strcmp('ps',src)
        % === Point source ===
        % Distance between loudspeaker and virtual source
        r = norm(x0-xs);
        % Delay and amplitude weight
        delay = r/c;
        weight = g0/(2*pi)*(x0-xs)*nx0'*r^(-3/2);

    elseif strcmp('fs',src)
        % === Focused source ===
        % Distance between loudspeaker and virtual source
        r = norm(x0-xs);
        % Delay and amplitude weight
        delay =  -r/c;
        weight = g0/(2*pi)*(x0-xs)*nx0'*r^(-3/2);
    else
        error('%s: %s is not a known source type.',upper(mfilename),src);
    end
else
    delay = 0;
    weight = 0;
end
