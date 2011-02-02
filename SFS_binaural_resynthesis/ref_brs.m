function brs = ref_brs(X,Y,phi,xs,ys,irs,conf)
%REF_BRS Generate a BRIR for a reference source
%   Usage: brs = ref_brs(X,Y,phi,xs,ys,irs,conf)
%          brs = ref_brs(X,Y,phi,xs,ys,irs)
%
%   Options:
%       X,Y     - listener position (m)
%       phi     - listener direction [head orientation] (°)
%       xs,ys   - source position (m)
%       irs     - IR data set for the second sources
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output:
%       brs     - Binaural room impulse response (nx2 matrix)
%
%   REF_BRS(X,Y,phi,xs,ys,irs,conf) calculates a binaural room impulse
%   response for a reference source (single loudspeaker) at position
%   [xs,ys] and a listener located at [X,Y].
%
%   Geometry:
%
%    x-axis
%       <---------------------------|----------------------------
%                                   |
%                 x [xs ys]         |
%           (Single Source)         |
%                                   |       |
%                                   |       O [X Y], phi
%                                   |    (Listener)
%                                   |
%                                   |
%                                   |
%                                   |
%                                   v y-axis
%
% see also: SFS_config, wfs_brs, auralize_brs, ref_brs_set
%

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input  parameters ==================================

if nargchk(6,7,nargin)
    error(['Wrong number of args.',...
           'Usage: brs = ref_brs(X,Y,phi,xs,ys,irs,conf)']);
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

if ~isstruct(irs)
    error('%s: irs has to be a struct!',upper(mfilename));
end

if nargin<7
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

Y0 = conf.Y0;                 % y coordinate of array

c = conf.c;                   % speed of sound

N = conf.N;                   % target length of BRS impulse responses

usehcomp = conf.usehcomp;     % Apply headphone compensation?
hcomplfile = conf.hcomplfile; % Headphone compensation file left
hcomprfile = conf.hcomprfile; % Headphone compensation file right


%% ===== variables ======================================================

% Check if we have a non focused source
if ys<Y0
    t0 = 0;
end

% Head orientation (counter clockwise, 0...2pi)
phi = phi/180*pi;

% HRIRs
lenir = length(irs.left(:,1));


%% ===== BRIR ===========================================================

% Initial values
brs = zeros(N,2);

% === Relative distances ===
%
% x-axis <---------------------------|---------------------------
%                                    |
%             [xs,ys]                |
%                x                   |
%                 |                  |
%                   | R              |
%                     |              |
%                       O            |
%                      [X,Y]         |
%                                    v
%                                  y-axis
%
% Distance between single loudspeaker [xs,ys] and listener position [X,Y]
R = norm( [X-xs, Y-ys] );

% Time delay of the single source (at the listener position)
% t0 is a causality pre delay for the focused source, e.g. 0 for a non
%>focused point source (see config.m)
% Single source
tau = ( (R-irs.distance)/c - t0 );
% Time delay in samples
dt = ceil( tau*fs );


% === Amplitude factor ===
% The 1/R term is for the decreasing of the sound on its way from
% the loudspeaker to the listener (R).
a = 1/R;


% === Secondary source angle ===
% Calculate the angle between the given loudspeaker and the listener.
% This is needed for the HRIR dataset.
%
% x-axis <---------------------------|---------------------------
%                       |            |
%             [xs,ys]   |            |
%                x   __ |            |
%                 | -   |            |  a = alpha
%                   | a |            |  cos(alpha) = dx(2)/R
%                 R   | |            |  tan(alpha) = -dx(1)/dx(2)
%                       O            |
%                      [X,Y]         |
%                                    v
%                                  y-axis
%
% Note: the above picture explains also that the cos(alpha) is not
% sufficient to span the whole number of possible angles! Only 0..180°
% is covered by acos!
% For the whole area tan2 is needed.
%
% Vector from listener position to single loudspeaker position (cp. R)
% NOTE: X - xs gives a negative value for a single loudspeaker to the left
% of the listener, therefor -dx1 is needed to get the right angle.
dx = [X Y] - [xs ys];
% Angle between listener and single loudspeaker source (-pi < alpha <= pi,
% without phi)
% Note: phi is the orientation of the listener (see first graph)
alpha = atan2(-dx(1),dx(2)) - phi;
%
% Ensure -pi <= alpha < pi
alpha = correct_azimuth(alpha);


% === HRIR interpolation ===
ir = get_ir(irs,alpha);

% Sum up virtual loudspeakers/HRIRs and add loudspeaker time delay
brs(:,1) = [zeros(1,dt) a*ir(:,1)' zeros(1,N-dt-lenir)]';
brs(:,2) = [zeros(1,dt) a*ir(:,2)' zeros(1,N-dt-lenir)]';


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

