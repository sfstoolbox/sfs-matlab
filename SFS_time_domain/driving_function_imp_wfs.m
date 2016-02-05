function [d,delay,weight] = driving_function_imp_wfs(x0,xs,src,conf)
%DRIVING_FUNCTION_IMP_WFS_25D calculates the WFS weighting and delaying
%
%   Usage: [d,delay,weight] = driving_function_imp_wfs(x0,xs,src,conf)
%
%   Input parameters:
%       x0      - positions and directions of secondary sources / m [nx6]
%       xs      - position of virtual source or direction of plane wave / m
%                 [1x3]
%       src     - source type of the virtual source
%                     'pw' - plane wave (xs, ys are the direction of the
%                            plane wave in this case)
%                     'ps' - point source
%                     'ls' - line source
%                     'fs' - focused source
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       d       - driving signals [mxn]
%       delay   - delay of the driving function / s [nx1]
%       weight  - weight (amplitude) of the driving function [nx1]
%
%   DRIVING_FUNCTION_IMP_WFS(x0,xs,src,conf) returns the driving signals and
%   weighting and delay parameters of the WFS driving function for the given
%   source type, position and secondary sources.
%
%   See also: sound_field_imp, sound_field_imp_wfs, driving_function_mono_wfs

%*****************************************************************************
% Copyright (c) 2010-2016 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2016 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
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
removedelay = conf.wfs.removedelay;


%% ===== Computation =====================================================
% Calculate pre-equalization filter if required
[pulse,prefilter_delay] = wfs_preequalization(dirac_imp(),conf);

% Secondary source positions and directions
nx0 = x0(:,4:6);
x0 = x0(:,1:3);

% Source position
xs = repmat(xs(1:3),[size(x0,1) 1]);

% Get the delay and weighting factors
if strcmp('pw',src)
    % === Plane wave =====================================================
    % Direction of plane wave
    nk = bsxfun(@rdivide,xs,vector_norm(xs,2));
    % Delay and amplitude weight
    [delay,weight] = driving_function_imp_wfs_pw(x0,nx0,nk,conf);

elseif strcmp('ps',src)
    % === Point source ===================================================
    % Delay and amplitude weight
    [delay,weight] = driving_function_imp_wfs_ps(x0,nx0,xs,conf);

elseif strcmp('ls',src)
    % === Line source ====================================================
    % Delay and amplitude weight
    [delay,weight] = driving_function_imp_wfs_ls(x0,nx0,xs,conf);

elseif strcmp('fs',src)
    % === Focused source =================================================
    % Delay and amplitude weight
    [delay,weight] = driving_function_imp_wfs_fs(x0,nx0,xs,conf);
else
    error('%s: %s is not a known source type.',upper(mfilename),src);
end

if removedelay
    % Remove delay offset, in order to begin always at t=0 with the first wave front
    % at any secondary source
    delay = delay-min(delay);
else
    % Delay to ensure causality at all secondary sources
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
    % Ensure virtual source start activity at t=0
    delay = delay + t0 + prefilter_delay;
end
% Append zeros at the end of the driving function. This is necessary, because
% the delayline function cuts into the end of the driving signals in order to
% delay them. NOTE: this can be changed by the conf.N setting
d_proto = repmat([row_vector(pulse) zeros(1,N-length(pulse))]',1,size(x0,1));
% Shift and weight prototype driving function
d = delayline(d_proto,delay*fs,weight,conf);
