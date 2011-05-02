function [weight,delay] = driving_function_imp_wfs_25d(x0,y0,phi,xs,ys,src,conf)
%DRIVING_FUNCTION_IMP_WFS_25D calculates the WFS 2.5D weighting and delaying
%   Usage: [weight,delay] = driving_function_imp_wfs_25d(...
%              x0,y0,phi,xs,ys,src,conf);
%          [weight,delay] = driving_function_imp_wfs_25d(x0,y0,phi,xs,ys,src);
%
%   Input parameters:
%       x0,y0   - position of secondary source (m)
%       phi     - orientation of secondary source (rad)
%       xs, ys  - position of virtual source or diirection of plane wave (m)
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
%   DRIVING_FUNCTION_IMP_WFS_25D(x0,y0,phi,xs,ys,src,conf) returns the
%   weighting and delay parameters of the WFS 2.5D driving function for the given
%   source type and position and loudspeaker positions.
%
%   see also: wave_field_imp_wfs_25d, driving_function_mono_wfs_25d

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 6;
nargmax = 7;
error(nargchk(nargmin,nargmax,nargin));
isargscalar(x0,y0,phi,xs,ys);
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
xref = conf.xref;
yref = conf.yref;


%% ===== Computation =====================================================

% Check also the activity of the used loudspeaker.
ls_activity = secondary_source_selection(x0,y0,phi,xs,ys,src);
% Constant amplitude factor
g0 = sqrt(2*pi*norm([xref yref]-[x0 y0]));
if ls_activity>0

    % Direction of secondary sources
    nx0 = -sin(phi);
    ny0 = cos(phi);

    if strcmp('pw',src)
        % === Plane wave ===

        % Direction of plane wave
        nxs = xs / sqrt(xs^2+ys^2);
        nys = ys / sqrt(xs^2+ys^2);
        theta = -1*atan2(nxs,nys);

        % delay in secs
        delay = (nxs*x0 + nys*y0)/c;
        % Amplitude
        % NOTE: <n_pw,n(x0)> is the same as the cosinus between their angle
        weight = 2*g0*cos(theta - phi);

    elseif strcmp('ps',src)
        % === Point source ===
        % Distance between loudspeaker and virtual source
        r = norm([x0 y0]-[xs ys]);
        % Delay and amplitude weight
        delay = r/c;
        weight = g0/(2*pi)*([x0 y0]-[xs ys])*[nx0 ny0]'*r^(-3/2);
    elseif strcmp('fs',src)
        % Focused source
        % Distance between loudspeaker and virtual source
        r = norm([x0 y0]-[xs ys]);
        % Delay and amplitude weight
        delay =  -r/c;
        weight = g0/(2*pi)*([x0 y0]-[xs ys])*[nx0 ny0]'*r^(-3/2);
    else
        error('%s: %s is not a known source type.',upper(mfilename),src);
    end
else
    delay = 0;
    weight = 0;
end
