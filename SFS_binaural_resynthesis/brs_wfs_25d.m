function brir = brs_wfs_25d(X,Y,phi,xs,ys,L,src,irs,conf)
%BRS_WFS_25D Generate a BRIR for WFS
%   Usage: brir = brs_wfs_25d(X,Y,phi,xs,ys,L,src,irs,conf)
%          brir = brs_wfs_25d(X,Y,phi,xs,ys,L,src,irs)
%
%   Input parameters:
%       X,Y     - listener position (m)
%       phi     - listener direction [head orientation] (rad)
%                 0 means the head is oriented towards the x-axis.
%       xs,ys   - virtual source position [ys > Y0 => focused source] (m)
%       L       - Length of linear loudspeaker array (m)
%       src     - source type: 'pw' -plane wave
%                              'ps' - point source
%                              'fs' - focused source
%       irs     - IR data set for the secondary sources
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameters:
%       brir    - Binaural room impulse response for the desired WFS array
%                 (nx2 matrix)
%
%   BRS_WFS_25D(X,Y,phi,xs,ys,L,irs,src,conf) calculates a binaural room impulse
%   response for a virtual source at [xs,ys] for a virtual WFS array and a
%   listener located at [X,Y].
%
%   Geometry:
%                               y-axis
%                                 ^
%                                 |
%                                 |
%                                 |    (Listener)
%                                 |        O [X Y], phi=-pi/2
%                                 |        |
%                                 |
%                  o [xs ys]      |
%             (Virtual Source)    |
%                                 |
%     -------v--v--v--v--v--v--v--v--v--v--v--v--v--v--v------> x-axis
%                              [X0 Y0] (Array center)
%            |---      Loudspeaker array length     ---|
%
% see also: brsset_wfs_25d, brs_point_source, auralize_brs
%

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 8;
nargmax = 9;
error(nargchk(nargmin,nargmax,nargin));

isargscalar(X,Y,phi,xs,ys);
isargpositivescalar(L);
isargchar(src);
check_irs(irs);

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================

fs = conf.fs;                 % sampling frequency
t0 = conf.t0;                 % pre-delay for causality (focused sources)
c = conf.c;                   % speed of sound
N = conf.N;                   % target length of BRS impulse responses
useplot = conf.useplot;       % Plot results?


%% ===== Variables ======================================================

phi = correct_azimuth(phi);

% Loudspeaker positions (phiLS describes the directions of the loudspeakers)
[x0,y0,phiLS] = secondary_source_positions(L,conf);
nls = length(x0);

% === Tapering window ===
win = tapwin(L,conf);

% === IRs ===
lenir = length(irs.left(:,1));


%% ===== BRIR ===========================================================

% Initial values
brir = zeros(N,2);
dt = zeros(1,nls);
a = zeros(1,nls);

% Create a BRIR for every single loudspeaker
warning('off','SFS:irs_intpol');
for n=1:nls

    % ====================================================================
    % Driving function to get weighting and delaying
    [a(n),delay] = ...
        driving_function_imp_wfs_25d(x0(n),y0(n),phiLS(n),xs,ys,src,conf);
    % Time delay of the virtual source (at the listener position)
    % t0 is a causality pre delay for a focused source, e.g. 0 for a non
    % focused point source (see SFS_config.m)
    % Define an offset to ensure |[X-Y]-[x0 y0]|-irs.distance+offset > 0
    offset = 3; % in m
    % Check if we have a non focused source
    if strcmp('fs',src)
        % Focused source
        tau = (norm([X Y]-[x0(n) y0(n)])-irs.distance+offset)/c + delay - t0;
    else
        % Virtual source behind the loudspeaker array
        tau = (norm([X Y]-[x0(n) y0(n)])-irs.distance+offset)/c + delay;
    end
    % Time delay in samples for the given loudspeaker
    % NOTE: I added some offset, because we can't get negative
    dt(n) = ceil( tau*fs ) + 300;
    if dt(n)<0
        error('%s: the time delay dt(n) = %i has to be positive.', ...
            upper(mfilename),dt(n));
    end
    if dt(n)<0
        error('%s: the time delay dt(n) = %i has to be positive.', ...
            upper(mfilename),dt(n));
    end

    % === Secondary source model: Greens function ===
    g = 1/(4*pi*norm([X Y]-[x0(n) y0(n)]));

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
    alpha = atan2(y0(n)-Y,x0(n)-X) - phi;
    %
    % Ensure -pi <= alpha < pi
    alpha = correct_azimuth(alpha);

    % === IR interpolation ===
    % Get the desired IR.
    % If needed interpolate the given IR set
    ir = get_ir(irs,alpha);

    % Check if we have enough samples (conf.N)
    if N<lenir+dt(n)
        error('Use a larger conf.N value, you need at least %i',lenir+dt(n));
    end

    % Sum up virtual loudspeakers/HRIRs and add loudspeaker time delay
    brir(:,1) = brir(:,1) + [zeros(1,dt(n)) ...
                             a(n)*win(n)*g*ir(:,1)' ...
                             zeros(1,N-dt(n)-lenir)]';
    brir(:,2) = brir(:,2) + [zeros(1,dt(n)) ...
                             a(n)*win(n)*g*ir(:,2)' ...
                             zeros(1,N-dt(n)-lenir)]';

    %figure;
    %title('IR');
    %plot([zeros(1,dt(n)) a(n)*win(n)*g*ir(:,1)' zeros(1,N-dt(n)-lenir)]');
    %figure;
    %title('BRIR');
    %plot(brir(:,1));

end
warning('on','SFS:irs_intpol');


%% ===== Pre-equalization ===============================================
brir = wfs_preequalization(brir,conf);

%% ===== Headphone compensation =========================================
brir = compensate_headphone(brir,conf);


%% ===== Plot WFS parameters ============================================
if(useplot)
    figure
    plot(x0,dt);
    title('delay (taps)');
    grid on;

    figure
    plot(x0,a);
    title('amplitude');
    grid on;
end
