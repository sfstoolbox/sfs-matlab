function brir = brs_point_source(X,Y,phi,xs,ys,irs,conf)
%BRS_POINT_SOURCE Generate a BRIR for a point source
%   Usage: brs = brs_point_source(X,Y,phi,xs,ys,irs,conf)
%          brs = brs_point_source(X,Y,phi,xs,ys,irs)
%
%   Options:
%       X,Y     - listener position (m)
%       phi     - listener direction [head orientation] (rad)
%       xs,ys   - source position (m)
%       irs     - IR data set for the second sources
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output:
%       brir    - Binaural room impulse response (nx2 matrix)
%
%   BRS_POINT_SOURCE(X,Y,phi,xs,ys,irs,conf) calculates a binaural room impulse
%   response for a reference source (single loudspeaker) at position
%   [xs,ys] and a listener located at [X,Y] and looking into direction phi.
%   Whereby at phi = 0 the listener is looking in the direction of the x-axis,
%   like the angle is normally defined in Mathematics.
%
%   Geometry:
%
%                                 y-axis
%                                   ^
%                                   |
%                                   |
%                                   |
%                                   |    listener
%                                   |       O [X Y], phi=-pi/2
%                                   |       |
%               source              |
%                 o [xs ys]         |
%                                   |
%                                   |
%       ----------------------------|---------------------------> x-axis
%
%
% see also: SFS_config, brs_wfs_25d, auralize_brs, brsset_point_source
%

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input parameters ====================================
nargmin = 6;
nargmax = 7;
error(nargchk(nargmin,nargmax,nargin));
isargscalar(X,Y,phi,xs,ys);
check_irs(irs);

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================

fs = conf.fs;                 % sampling frequency
t0 = conf.t0;                 % pre-delay for causality (focused sources)
Y0 = conf.Y0;                 % y coordinate of array
c = conf.c;                   % speed of sound
N = conf.N;                   % target length of BRS impulse responses
usehcomp = conf.usehcomp;     % Apply headphone compensation?
hcomplfile = conf.hcomplfile; % Headphone compensation file left
hcomprfile = conf.hcomprfile; % Headphone compensation file right


%% ===== variables =======================================================

phi = correct_azimuth(phi);

% HRIRs
lenir = length(irs.left(:,1));


%% ===== BRIR ============================================================

% Initial values
brir = zeros(N,2);

% === Secondary source angle ===
% Calculate the angle between the given loudspeaker and the listener.
% This is needed for the HRIR dataset.
%
%                                 y-axis
%                                    ^
%                                    |
%               [X,Y], phi = 0       |
%                 O--------          |
%                  \ a |             |
%                   \ /              |  a = alpha
%                    \               |  tan(alpha) = (ys-Y)/(xs-X)
%                     o              |
%                  [xs,ys]           |
%        ----------------------------|--------------------------> x-axis
%
% Angle between listener and source (-pi < alpha <= pi)
% NOTE: phi is the orientation of the listener (see first graph)
alpha = atan2(ys-Y,xs-X) - phi;
%
% Ensure -pi <= alpha < pi
alpha = correct_azimuth(alpha);

% === HRIR interpolation ===
ir = get_ir(irs,alpha);
ir_distance = get_ir_distance(irs,alpha);

% === Relative distances and time delay ===
% Distance between single loudspeaker [xs,ys] and listener position [X,Y]
R = norm([X Y]-[xs ys]);
% Time delay of the single source (at the listener position)
% Define an offset to ensure R-irs.distance+offset > 0
offset = 3; % in m
tau = (R-ir_distance+offset)/c;
% Time delay in samples
dt = ceil( tau*fs );

% === Amplitude factor ===
% The 1/R term is for the decreasing of the sound on its way from
% the loudspeaker to the listener (R). It accounts for the distance that is
% already present in the IR dataset.
a = (1/R) / (1/ir_distance);

% Check if we have enough samples (conf.N)
if N<lenir+dt
    error('Use a larger conf.N value, you need at least %i',lenir+dt);
end
% Sum up virtual loudspeakers/HRIRs and add loudspeaker time delay
brir(:,1) = [zeros(1,dt) a*ir(:,1)' zeros(1,N-dt-lenir)]';
brir(:,2) = [zeros(1,dt) a*ir(:,2)' zeros(1,N-dt-lenir)]';


%% ===== Headphone compensation ==========================================
brir = compensate_headphone(brir,conf);
