function brs = wfs_brs(X,Y,phi,xs,ys,L,irs,conf)
%WFS_BRS Generate a BRIR for WFS
%   Usage: brs = wfs_brs(X,Y,phi,xs,ys,L,irs,conf)
%          brs = wfs_brs(X,Y,phi,xs,ys,L,irs)
%
%   Input parameters:
%       X,Y     - listener position (m)
%       phi     - listener direction [head orientation] (�)
%       xs,ys   - virtual source position [ys > Y0 => focused source] (m)
%       L       - Length of linear loudspeaker array (m)
%       irs     - IR data set for the secondary sources
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameters:
%       brs     - Binaural room impulse response (nx2 matrix)
%
%   WFS_BRS(X,Y,phi,xs,ys,L,irs,conf) calculates a binaural room impulse
%   response for a virtual source at [xs,ys] for a linear WFS array and the
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
% see also: SFS_config, LSpos_linear, ref_brs, auralize_brs, wfs_brs_set
%

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input  parameters ==================================

if nargchk(7,8,nargin)
    error(['Wrong number of args.',...
           'Usage: brs = wfs_brs(X,Y,phi,xs,ys,L,irs,conf)']);
end

if ~isnumeric(X) || ~isscalar(X)
    error('%s: X has to be a scalar!',upper(mfilename));
end

if ~isnumeric(Y) || ~isscalar(Y)
    error('%s: Y has to be a scalar!',upper(mfilename));
end

if ~isnumeric(phi) || ~isscalar(phi)
    error('%s: phi has to be a scalar!',upper(mfilename));
end

if ~isnumeric(xs) || ~isscalar(xs)
    error('%s: xs has to be a scalar!',upper(mfilename));
end

if ~isnumeric(ys) || ~isscalar(ys)
    error('%s: ys has to be a scalar!',upper(mfilename));
end

if ~isnumeric(L) || ~isscalar(L) || L<0
    error('%s: L has to be a positive scalar!',upper(mfilename));
end

if ~isstruct(irs)
    error('%s: irs has to be a struct!',upper(mfilename));
end

if nargin<8
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

fs = conf.fs;                 % sampling frequency
t0 = conf.t0;                 % pre-delay for causality (focused sources)

LSdist = conf.LSdist;         % loudspeaker distance
c = conf.c;                   % speed of sound

X0 = conf.X0;                 % array position
Y0 = conf.Y0;

N = conf.N;                   % target length of BRS impulse responses

usehcomp = conf.usehcomp;     % Apply headphone compensation?
hcomplfile = conf.hcomplfile; % Headphone compensation file left
hcomprfile = conf.hcomprfile; % Headphone compensation file right

usehpre = conf.usehpre;       % Apply WFS pre-filter?

useplot = conf.useplot;       % Plot results?

% Check if the listener position is in front of the loudspeaker array
if Y<=Y0
    error(['%s: The listener position Y is located behind the'...
           'loudspeaker (Y0)'],upper(mfilename));
end


%% ===== Variables ======================================================

% Head orientation (counter clockwise, 0...2pi)
phi = phi/180*pi;

% Number of loudspeakers (round towards plus infinity)
nLS = ceil(L/LSdist);

% Loudspeaker positions (LSdir describes the directions of the LS) for a
% linear WFS array
[LSpos,LSdir] = LSpos_linear(X0,Y0,(nLS-1)*LSdist,nLS);

% === Tapering window ===
% See in config.m if it is applied
win = tapwin(L,conf);

% === HRIRs ===
lenir = length(irs.left(:,1));


%% ===== BRIR ===========================================================

% Initial values
brs = zeros(N,2);
dt = zeros(1,nLS);
a = zeros(1,nLS);

% Create a BRIR for every single loudspeaker
for n=1:nLS

    % === Relative distances ===
    %
    %           LSpos(:,n)            [X0 Y0]
    % x-axis <-^--^--^--^--^--^--^--^--^-|-^--^--^--^--^--^--^--^--^--
    %             |                      |
    %         R2 |  |                    |
    %           |     | R                |
    %          x        |                |
    %       [xs,ys]       |              |
    %                       O            |
    %                      [X,Y]         |
    %                                    v
    %                                  y-axis

    % Distance between given loudspeaker and listener position [X,Y]
    R = norm( [X-LSpos(1,n), Y-LSpos(2,n)] );
    % Distance between given loudspeaker and virtual source position
    R2 = norm( [xs-LSpos(1,n), ys-LSpos(2,n)] );

    % Time delay of the virtual source (at the listener position)
    % r0 is the measurement distance for the HRIR data (see SFS_config.m)
    % t0 is a causality pre delay for the focused source, e.g. 0 for a non
    % focused point source (see SFS_config.m)
    % Check if we have a non focused source
    if ys>Y0
        % Focused source
        tau = (R-irs.r0)/c - R2/c - t0;
    else
        % Virtual source behind the loudspeaker array
        tau = (R-irs.r0)/c + R2/c;
    end
    % Time delay in samples for the given loudspeaker
    dt(n) = ceil( tau*fs );


    % === Amplitude factor ===
    % Use linear amplitude factor (see Spors et al. (2008)) and apply the
    % tapering window
    a(n) = wfs_amplitude_linear(LSpos(1,n),LSpos(2,n),X,Y,xs,ys) * win(n);


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
    % sufficient to span the whole number of possible angles! Only 0..180�
    % is covered by acos!
    % For the whole area tan2 is needed.
    %
    % Vector from listener position to given loudspeaker position (cp. R)
    % NOTE: X - LSpos(1,n) gives a negative value for loudspeaker to the
    % left of the listener, therefor -dx1 is needed to get the right angle.
    dx = [X Y]' - LSpos(:,n);
    % Angle between listener and secondary source (-pi < alpha <= pi,
    % without phi)
    % Note: phi is the orientation of the listener (see first graph)
    alpha = atan2(-dx(1),dx(2)) - phi;
    %
    % Ensure -pi <= alpha < pi
    alpha = correct_angle(alpha);

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
    brs(:,1) = brs(:,1) + [zeros(1,dt(n)) ...
                             a(n)*ir(:,1)' ...
                             zeros(1,N-dt(n)-lenir)]';
    brs(:,2) = brs(:,2) + [zeros(1,dt(n)) ...
                             a(n)*ir(:,2)' ...
                             zeros(1,N-dt(n)-lenir)]';

end


%% ===== Pre-equalization ===============================================
if(usehpre)
    % Generate WFS preequalization-filter
    hpre = wfs_prefilter(conf);
    % Apply filter
    brs(:,1) = conv(hpre,brs(1:end-length(hpre)+1,1));
    brs(:,2) = conv(hpre,brs(1:end-length(hpre)+1,2));
end


%% ===== Headphone compensation =========================================
if(usehcomp)
    % Read headphone compensation filter
    hcompl = wavread(hcomplfile);
    hcompr = wavread(hcomprfile);
    hcomp = [hcompl hcompr];
    % Apply filter
    brs(:,1) = conv(hcomp(:,1),brs(1:end-length(hcomp)+1,1));
    brs(:,2) = conv(hcomp(:,2),brs(1:end-length(hcomp)+1,2));
end


%% ===== Plot WFS parameters ============================================
if(useplot)
    figure
    plot(LSpos(1,:),dt);
    title('delay (taps)');
    grid on;

    figure
    plot(LSpos(1,:),a);
    title('amplitude');
    grid on;
end

