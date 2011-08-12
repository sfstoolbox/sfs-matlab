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
%              |---      Loudspeaker array length     ---|
%    x-axis                      [X0 Y0] (Array center)
%       <------^--^--^--^--^--^--^--^--^--^--^--^--^--^--^-------
%                                   |
%                 x [xs ys]         |
%           (Virtual Source)        |
%                                   |        |
%                                   |        O [X Y], phi
%                                   |    (Listener)
%                                   |
%                                   |
%                                   |
%                                   |
%                                   v y-axis
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

% Loudspeaker positions (LSdir describes the directions of the LS) for a
% linear WFS array
[x0,y0,phiLS] = secondary_source_positions(L,conf);
nLS = length(x0);


% === Tapering window ===
% See in SFS_config.m if it is applied
win = tapwin(L,conf);

% === HRIRs ===
lenir = length(irs.left(:,1));


%% ===== BRIR ===========================================================

% Initial values
brir = zeros(N,2);
dt = zeros(1,nLS);
a = zeros(1,nLS);

% Create a BRIR for every single loudspeaker
warning('off','SFS:irs_intpol');
for n=1:nLS

    % ====================================================================
    % Driving function to get weighting and delaying
    [a(n),delay] = ...
        driving_function_imp_wfs_25d(x0(n),y0(n),phiLS(n),xs,ys,src,conf);
    % Time delay of the virtual source (at the listener position)
    % t0 is a causality pre delay for the focused source, e.g. 0 for a non
    % focused point source (see SFS_config.m)
    % Check if we have a non focused source
    if strcmp('fs',src)
        % Focused source
        %tau = (norm([X Y]-[x0(n) y0(n)]) - irs.distance)/c + delay - t0;
        tau = norm([X Y]-[x0(n) y0(n)])/c + delay - t0;
    else
        % Virtual source behind the loudspeaker array
        %tau = (norm([X Y]-[x0(n) y0(n)])-irs.distance)/c + delay;
        tau = norm([X Y]-[x0(n) y0(n)])/c + delay;
    end
    % Time delay in samples for the given loudspeaker
    % NOTE: I added some offset, because we can't get negative
    dt(n) = ceil( tau*fs ) + 300;
    if dt(n)<0
        error('%s: the time delay dt(n) = %i has to be positive.', ...
            upper(mfilename),dt(n));
    end


    % === Secondary source model: Greens function ===
    g = 1/(4*pi*norm([x0(n) y0(n)]-[X Y]));


    % === Secondary source angle ===
    % Calculate the angle between the given loudspeaker and the listener.
    % This is needed for the HRIR dataset.
    %
    %          LSpos -dx(1)           [X0 Y0]
    % x-axis <-^--^--^--^--^--^--^--^--^-|-^--^--^--^--^--^--^--^--^-
    %             |         |            |
    %               |       |            |
    %               R | _---| dx(2)      |
    %                   | a |            |  a = alpha
    %                     | |            |  cos(alpha) = dx(2)/R
    %                       O            |  tan(alpha) = -dx(1)/dx(2)
    %                      [X,Y]         |
    %                                    v
    %                                  y-axis
    %
    % Note: the above picture explains also that the cos(alpha) is not
    % sufficient to span the whole number of possible angles! Only 0..180 deg
    % is covered by acos!
    % For the whole area tan2 is needed.
    %
    % Vector from listener position to given loudspeaker position (cp. R)
    % NOTE: X - x0(n) gives a negative value for loudspeaker to the
    % left of the listener, therefor -dx1 is needed to get the right angle.
    dx = [X Y] - [x0(n) y0(n)];
    % Angle between listener and secondary source (-pi < alpha <= pi,
    % without phi)
    % Note: phi is the orientation of the listener (see first graph)
    alpha = atan2(-dx(1),dx(2)) - phi;
    %
    % Ensure -pi <= alpha < pi
    alpha = correct_azimuth(alpha);

    % === IR interpolation ===
    % Get the desired IR.
    % If needed interpolate the given IR set
    ir = get_ir(irs,alpha);

    % Check if the length of the BRIR (conf.N) is long enough for the
    % needed time delay dt(n)
    if N<(dt(n)+lenir)
        error(['%s: The length of the BRIR conf.N is not large enough ' ...
               'to handle the needed time delay dt(n).'],upper(mfilename));
    end

    % Sum up virtual loudspeakers/HRIRs and add loudspeaker time delay
    brir(:,1) = brir(:,1) + [zeros(1,dt(n)) ...
                             a(n)*win(n)*g*ir(:,1)' ...
                             zeros(1,N-dt(n)-lenir)]';
    brir(:,2) = brir(:,2) + [zeros(1,dt(n)) ...
                             a(n)*win(n)*g*ir(:,2)' ...
                             zeros(1,N-dt(n)-lenir)]';

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

