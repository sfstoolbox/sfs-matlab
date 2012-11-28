function ir_wfs = ir_wfs_25d(X,phi,xs,src,L,irs,conf)
%BRS_WFS_25D Generate a IR for WFS
%
%   Usage: ir_wfs = ir_wfs_25d(X,phi,xs,src,L,irs,[conf])
%
%   Input parameters:
%       X       - listener position (m)
%       phi     - listener direction [head orientation] (rad)
%                 0 means the head is oriented towards the x-axis.
%       xs      - virtual source position [ys > Y0 => focused source] (m)
%       src     - source type: 'pw' -plane wave
%                              'ps' - point source
%                              'fs' - focused source
%       L       - Length of loudspeaker array (m)
%       irs     - IR data set for the secondary sources
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       iri_wfs - Impulse response for the desired WFS array (nx2 matrix)
%
%   IR_WFS_25D(X,phi,xs,src,L,irs,conf) calculates a binaural room impulse
%   response for a virtual source at xs for a virtual WFS array and a
%   listener located at X.
%
%   Geometry (for a linear array):
%
%                               y-axis
%                                 ^
%                                 |
%                                 |
%                                 |    (Listener)
%                                 |        O X, phi=-pi/2
%                                 |        |
%                                 |
%                  o xs           |
%             (Virtual Source)    |
%                                 |
%     -------v--v--v--v--v--v--v--v--v--v--v--v--v--v--v------> x-axis
%                                 X0 (Array center)
%            |---      Loudspeaker array length     ---|
%
% see also: brs_wfs_25d, ir_point_source, auralize_ir

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 6;
nargmax = 7;
narginchk(nargmin,nargmax);
if nargin==nargmax-1
    conf = SFS_config;
end
[X,xs] = position_vector(X,xs);
if conf.debug
    isargscalar(phi);
    isargpositivescalar(L);
    isargchar(src);
    check_irs(irs);
end


%% ===== Configuration ==================================================
fs = conf.fs;                 % sampling frequency
t0 = conf.t0;                 % pre-delay for causality (focused sources)
c = conf.c;                   % speed of sound
N = conf.N;                   % target length of BRS impulse responses
xref = conf.xref;             % reference point (needed for direction of focused
                              % sources)
debug = conf.debug;           % debugging


%% ===== Variables ======================================================

phi = correct_azimuth(phi);

% Loudspeaker positions (phiLS describes the directions of the loudspeakers)
x0 = secondary_source_positions(L,conf);
x0 = secondary_source_selection(x0,xs,src,xref);
nls = size(x0,1);
% generate tapering window
win = tapering_window(x0,conf);


%% ===== BRIR ===========================================================
% Initial values
ir_wfs = zeros(N,2);
dt = zeros(1,nls);
a = zeros(1,nls);

% Create a BRIR for every single loudspeaker
warning('off','SFS:irs_intpol');
for ii = 1:nls


    % === Secondary source model: Greens function ===
    g = 1/(4*pi*norm(X-x0(ii,1:3)));

    % === Secondary source angle ===
    % Calculate the angle between the given loudspeaker and the listener.
    % This is needed for the HRIR dataset.
    %
    %                             y-axis
    %                               ^
    %                               |
    %                               |
    %                               |
    %            [X Y], phi=0       |
    %              O------------    |  a = alpha
    %               \ a |           |  tan(alpha) = (y0-Y)/(x0-X)
    %                \ /            |
    %                 \             |
    %                  \            |
    %   -------v--v--v--v--v--v--v--v--v--v--v--v--v--v--v------> x-axis
    %                [x0 y0]
    %
    % Angle between listener and secondary source (-pi < alpha <= pi)
    % Note: phi is the orientation of the listener (see first graph)
    [alpha,theta_tmp,r_tmp] = cart2sph(x0(ii,1)-X(1),x0(ii,2)-X(2),0);
    %
    % Ensure -pi <= alpha < pi
    alpha = correct_azimuth(alpha-phi);

    % === IR interpolation ===
    % Get the desired IR.
    % If needed interpolate the given IR set
    ir = get_ir(irs,alpha,0);
    ir_distance = get_ir_distance(irs,alpha,0);

    % === Amplitude and delay ===
    % Driving function to get weighting and delaying
    [a(ii),delay] = driving_function_imp_wfs_25d(x0(ii,:),xs,src,conf);
    % Time delay of the virtual source (at the listener position)
    % t0 is a causality pre delay for a focused source, e.g. 0 for a non
    % focused point source (see SFS_config.m)
    % Define an offset to ensure |[X-Y]-[x0 y0]|-irs.distance+offset > 0
    offset = 3; % in m
    % Check if we have a non focused source
    % FIXME: added a hack to handle the time gap in plane waves (this occurs
    % only for some directions). Find a better solution for this
    if strcmp('fs',src) || strcmp('pw',src)
        % Focused source
        tau = (norm(X-x0(ii,1:3)) - ir_distance+offset)/c + delay - t0;
    else
        % Virtual source behind the loudspeaker array
        tau = (norm(X-x0(ii,1:3)) - ir_distance+offset)/c + delay;
    end
    % Time delay in samples for the given loudspeaker
    % NOTE: I added some offset, because we can't get negative
    dt(ii) = (tau*fs) + 0;
    if dt(ii)<0
        error('%s: the time delay dt(ii) = %i has to be positive.', ...
            upper(mfilename),dt(ii));
    end

    % === Trim IR ===
    ir = fix_ir_length(ir,N,dt(ii));

    % === Sum up virtual loudspeakers/HRIRs and add loudspeaker time delay ===
    ir_wfs(:,1) = ir_wfs(:,1) + delayline(ir(:,1)',dt(ii),a(ii)*win(ii)*g,conf)';
    ir_wfs(:,2) = ir_wfs(:,2) + delayline(ir(:,2)',dt(ii),a(ii)*win(ii)*g,conf)';

end
warning('on','SFS:irs_intpol');


%% ===== Pre-equalization ===============================================
ir_wfs = wfs_preequalization(ir_wfs,conf);

%% ===== Headphone compensation =========================================
ir = compensate_headphone(ir_wfs,conf);


%% ===== Plot WFS parameters ============================================
if(debug)
    figure
    plot(x0(:,1),dt);
    title('delay (taps)');
    grid on;

    figure
    plot(x0(:,1),a);
    title('amplitude');
    grid on;
end
