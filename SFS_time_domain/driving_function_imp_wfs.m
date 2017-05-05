function [d,delay,weight,delay_offset] = driving_function_imp_wfs(x0,xs,src,conf)
%DRIVING_FUNCTION_IMP_WFS_25D calculates the WFS weighting and delaying
%
%   Usage: [d,delay,weight] = driving_function_imp_wfs(x0,xs,src,conf)
%
%   Input parameters:
%       x0              - positions and directions of secondary
%                         sources / m [nx7]
%       xs              - position of virtual source or direction of
%                         plane wave / m [1x3]
%       src             - source type of the virtual source
%                           'pw' - plane wave (xs, ys are the direction of the
%                                  plane wave in this case)
%                           'ps' - point source
%                           'ls' - line source
%                           'fs' - focused source
%       conf            - configuration struct (see SFS_config)
%
%   Output parameters:
%       d               - driving signals [mxn]
%       delay           - delay of the driving function / s [nx1]
%       weight          - weight (amplitude) of the driving function [nx1]
%       delay_offset    - additional added delay / s
%
%   DRIVING_FUNCTION_IMP_WFS(x0,xs,src,conf) returns the driving signals and
%   weighting and delay parameters of the WFS driving function for the given
%   source type, position and secondary sources.
%
%   See also: sound_field_imp, sound_field_imp_wfs, driving_function_mono_wfs

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargsecondarysource(x0)
isargxs(xs);
isargchar(src);
isargstruct(conf);


%% ===== Configuration ==================================================
fs = conf.fs;
N = conf.N;
c = conf.c;
t0 = conf.t0;


%% ===== Computation =====================================================
% Calculate pre-equalization filter if required
[pulse,prefilter_delay] = wfs_preequalization(dirac_imp(),conf);

% Secondary source positions and directions
nx0 = x0(:,4:6);
x0 = x0(:,1:3);

% Source position
xs = repmat(xs,[size(x0,1) 1]);

% Get the delay and weighting factors
if strcmp('pw',src)
    % === Plane wave ===
    % Direction of plane wave
    nk = bsxfun(@rdivide,xs,vector_norm(xs(:,1:3),2));
    % Delay and amplitude weight
    [delay,weight] = driving_function_imp_wfs_pw(x0,nx0,nk,conf);

elseif strcmp('ps',src)
    % === Point source ===
    [delay,weight] = driving_function_imp_wfs_ps(x0,nx0,xs(:,1:3),conf);

elseif strcmp('ls',src)
    % === Line source ===
    [delay,weight] = driving_function_imp_wfs_ls(x0,nx0,xs,conf);

elseif strcmp('fs',src)
    % === Focused source ===
    [delay,weight] = driving_function_imp_wfs_fs(x0,nx0,xs(:,1:3),conf);
else
    error('%s: %s is not a known source type.',upper(mfilename),src);
end

if strcmp('system',t0)
    % Set minimum delay to 0, in order to begin always at t=1 with the first
    % wave front at any secondary source
    delay = delay - min(delay);
    % Return extra offset due to prefilter
    delay_offset = prefilter_delay;
elseif strcmp('source',t0)
    % Add extra delay to ensure causality at all secondary sources (delay>0)
    [diameter,center] = secondary_source_diameter(conf);
    t0 = diameter/c;
    if (ceil((max(delay)+t0)*fs) - 1 ) > N
        % This is a lot more likely to happen.
        warning('conf.N = %i is too short for requested source.',N);
    end
    if strcmp('fs',src)
        % Reject focused source that's too far away
        % (this will only happen for unbounded listening arrays.)
        if norm(xs(1:3) - center,2) > diameter/2;
            error(['%s: Using ''config.wfs.removedelay == 0'', ', ...
                'focused source positions are restricted to the ball ', ...
                'spanned by the array diameter.'],upper(mfilename));
        end
    end
    delay = delay + t0;
    % Return extra added delay. This is can be used to ensure that the virtual
    % source always starts at t = 0.
    delay_offset = t0 + prefilter_delay;
else
    error('%s: t0 needs to be "system" or "source" and not "%s".', ...
          upper(mfilename),t0);
end
% Append zeros at the end of the driving function. This is necessary, because
% the delayline function cuts into the end of the driving signals in order to
% delay them. NOTE: this can be changed by the conf.N setting
d_proto = repmat([row_vector(pulse) zeros(1,N-length(pulse))]',1,size(x0,1));
% Shift and weight prototype driving function
[d, delayline_delay] = delayline(d_proto,delay,weight,conf);
% Add delay offset of delayline
delay_offset = delay_offset + delayline_delay;
