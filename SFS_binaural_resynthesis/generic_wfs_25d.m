function ir = generic_wfs_25d(xs,L,src,conf)
%GENRIC_WFS_25D Generate a IR for the generic renderer of the SSR
%
%   Usage: ir = generic_wfs_25d(xs,L,src,conf)
%          ir = generic_wfs_25d(xs,L,src)
%
%   Input parameters:
%       xs      - virtual source position [ys > Y0 => focused source] (m)
%       L       - Length of linear loudspeaker array (m)
%       src     - source type: 'pw' -plane wave
%                              'ps' - point source
%                              'fs' - focused source
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameters:
%       ir      - Impulse response for the desired WFS array loudspeaker
%                 (nx1)
%
%   GENERIC_WFS_25D(xs,L,src,conf) calculates an impulse
%   response for a virtual source at xs for the loudspeakers of a WFS
%   array. every loudspeaker of the array is represented by one column in
%   the impulse response.
%
% see also: brs_wfs_25d, brs_point_source, auralize_ir
%

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 5;
error(nargchk(nargmin,nargmax,nargin));
xs = position_vector(xs);
isargpositivescalar(L);
isargchar(src);

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
fs = conf.fs;                 % sampling frequency
c = conf.c;                   % speed of sound
N = conf.N;                   % target length of BRS impulse responses
useplot = conf.useplot;       % Plot results?


%% ===== Variables ======================================================
% Secondary sources
x0 = secondary_source_positions(L,conf);
nls = length(x0);
ls_activity = secondary_source_selection(x0,xs,src);
% generate tapering window
win = tapwin(L,ls_activity,conf);


%% ===== IR =============================================================

% === Dirac pulse ===
% Generate a dirac pulse to use in the delaying and add procedure
dirac = zeros(1024,1);
dirac(300) = 1;
lenir = length(dirac);

% Initial values
ir = zeros(N,nls);
dt = zeros(1,nls);
a = zeros(1,nls);

% Create a IR for every single loudspeaker
for n=1:nls

    % ====================================================================
    % Driving function to get weighting and delaying
    [a(n),delay] = ...
        driving_function_imp_wfs_25d(x0(n,:),xs,src,conf);
    % Time delay in samples for the given loudspeaker
    dt(n) = ceil( delay*fs ) + 500;

    % Check if the length of the IR (conf.N) is long enough for the
    % needed time delay dt(n)
    if N<(dt(n)+lenir)
        error(['%s: The length of the IR conf.N is not large enough ' ...
               'to handle the needed time delay dt(n).'],upper(mfilename));
    end

    % Generate impulse response for the desired loudspeaker
    ir(:,n) = [zeros(1,dt(n)) a(n)*win(n)*dirac' zeros(1,N-dt(n)-lenir)]';

end

%% ===== Pre-equalization ================================================
ir = wfs_preequalization(ir,conf);


%% ===== Plot WFS parameters =============================================
if(useplot)
    figure
    plot(x0(:,1),dt);
    title('delay (taps)');
    grid on;

    figure
    plot(x0(:,1),a);
    title('amplitude');
    grid on;
end

