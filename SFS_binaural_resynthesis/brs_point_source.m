function brir = brs_point_source(X,phi,xs,irs,conf)
%BRS_POINT_SOURCE Generate a BRIR for a point source
%   Usage: brs = brs_point_source(X,phi,xs,irs,conf)
%          brs = brs_point_source(X,phi,xs,irs)
%
%   Options:
%       X       - listener position (m)
%       phi     - listener direction [head orientation] (rad)
%       xs      - source position (m)
%       irs     - IR data set for the second sources
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output:
%       brir    - Binaural room impulse response (nx2 matrix)
%
%   BRS_POINT_SOURCE(X,phi,xs,irs,conf) calculates a binaural room impulse
%   response for a reference source (single loudspeaker) at position
%   xs and a listener located at X and looking into direction phi.
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
%                                   |       O X, phi=-pi/2
%                                   |       |
%               source              |
%                 o xs              |
%                                   |
%                                   |
%       ----------------------------|---------------------------> x-axis
%
%
% see also: SFS_config, brs_wfs_25d, auralize_ir, brsset_point_source
%

% AUTHOR: Sascha Spors, Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input parameters ====================================
nargmin = 4;
nargmax = 5;
error(nargchk(nargmin,nargmax,nargin));
[X,xs] = position_vector(X,xs);
isargscalar(phi);
check_irs(irs);

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
fs = conf.fs;                  % sampling frequency
t0 = conf.t0;                  % pre-delay for causality (focused sources)
X0 = position_vector(conf.X0); % center coordinate of array
c = conf.c;                    % speed of sound
N = conf.N;                    % target length of BRS impulse responses
usehcomp = conf.usehcomp;      % Apply headphone compensation?
hcomplfile = conf.hcomplfile;  % Headphone compensation file left
hcomprfile = conf.hcomprfile;  % Headphone compensation file right


%% ===== Variables =======================================================
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
%                X, phi = 0          |
%                 O--------          |
%                  \ a |             |
%                   \ /              |  a = alpha
%                    \               |
%                     o              |
%                     xs             |
%        ----------------------------|--------------------------> x-axis
%
% Angle between listener and source (-pi < alpha <= pi)
% NOTE: phi is the orientation of the listener (see first graph)
alpha = cart2sph(xs(1)-X(1),xs(2)-X(2),0) - phi;
%
% Ensure -pi <= alpha < pi
alpha = correct_azimuth(alpha);

% === HRIR interpolation ===
ir = get_ir(irs,alpha);
ir_distance = get_ir_distance(irs,alpha);

% === Relative distances and time delay ===
% Time delay of the single source (at the listener position)
% Define an offset to ensure norm(X-xs)-irs.distance+offset > 0
offset = 3; % in m
tau = (norm(X-xs)-ir_distance+offset)/c;
% Time delay in samples
dt = ceil( tau*fs );

% === Amplitude factor ===
% The 1/norm(X-xs) term is for the decreasing of the sound on its way from
% the loudspeaker (xs) to the listener (X). It accounts for the distance that is
% already present in the IR dataset.
a = (1/norm(X-xs)) / (1/ir_distance);

% FIXME: this could go in an extra function
% append zeros or truncate IRs to target length
if(lenir<N-dt)
    ir=cat(1,ir,zeros(N-lenir,2));
else
    ir=ir(1:N-dt,:);
end
% Check if we have enough samples (conf.N)
if N<lenir+dt
    %error('Use a larger conf.N value, you need at least %i',lenir+dt);
end
% Sum up virtual loudspeakers/HRIRs and add loudspeaker time delay
brir(:,1) = [zeros(1,dt) a*ir(:,1)' zeros(1,N-dt-lenir)]';
brir(:,2) = [zeros(1,dt) a*ir(:,2)' zeros(1,N-dt-lenir)]';


%% ===== Headphone compensation ==========================================
brir = compensate_headphone(brir,conf);
