function [f,S] = freq_response_wfs_25d(X,xs,L,src,conf)
%FREQ_RESPONSE_WFS_25D simulates the frequency response for 2.5D WFS
%
%   Usage: [f,S] = freq_response_wfs_25d(X,xs,L,src,conf)
%          [f,S] = freq_response_wfs_25d(X,xs,L,src)
%
%   Input parameters:
%       X           - listener position (m)
%       xs          - position of virtual source (m)
%       L           - array length (m)
%       src         - source type of the virtual source
%                         'pw' -plane wave
%                         'ps' - point source
%                         'fs' - focused source
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       f           - corresponding frequency (x) axis
%       S           - simulated frequency response
%
%   FREQ_RESPONSE_WFS_25D(X,xs,L,src,conf) simulates the frequency 
%   response of the wave field at the given position X. The wave field is 
%   simulated for the given source type (src) using a WFS 2.5 dimensional 
%   driving function in the temporal domain.
%
%   References:
%       Spors2009 - Physical and Perceptual Properties of Focused Sources in
%           Wave Field Synthesis (AES127)
%       Spors2010 - Analysis and Improvement of Pre-equalization in
%           2.5-Dimensional Wave Field Synthesis (AES128)
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: wave_field_mono_wfs_25d, wave_field_time_domain_wfs_25d
%

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 5;
error(nargchk(nargmin,nargmax,nargin));
[X,xs] = position_vector(X,xs);
isargpositivescalar(L);
isargchar(src);

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
% Array position (m)
X0 = position_vector(conf.X0);
% Reference position for the amplitude (correct reproduction of amplitude
% at y = yref).
yref = conf.yref;
% Use tapering window?
usetapwin = conf.usetapwin;
% xy resolution
xysamples = conf.xysamples;
% Phase of the wave field
phase = conf.phase;
% Plotting result
useplot = conf.useplot;
% Speed of sound
c = conf.c;


%% ===== Computation ====================================================
% Calculate the wave field in time domain

% Get the position of the loudspeakers
x0 = secondary_source_positions(L,conf);
nls = size(x0,1);

% Generate frequencies (10^0-10^5)
f = logspace(0,5,500);
% We want only frequencies until f = 20000Hz
idx = find(f>20000,1);
f = f(1:idx);

S = zeros(1,length(f));

% Activity of secondary sources
ls_activity = secondary_source_selection(x0,xs,src);
% Tapering window
win = tapwin(L,ls_activity,conf);
% Get the result for all frequencies
for ii = 1:length(f)
    P = 0;
    % Integration over secondary source positions
    for n = 1:nls

        % ================================================================
        % Secondary source model
        % This is the model for the loudspeakers we apply. We use closed cabinet
        % loudspeakers and therefore the 3D Green's function is our model.
        G = point_source(X(1),X(2),x0(n,1:3),f(ii),conf);

        % ================================================================
        % Driving function D(x0,omega)
        D = driving_function_mono_wfs_25d(x0(n,:),xs,f(ii),src,conf);

        % ================================================================
        % Integration
        %              /
        % P(x,omega) = | D(x0,omega) G(x-xs,omega) dx0
        %              /
        %
        % see: Spors2009, Williams1993 p. 36
        %
        P = P + win(n)*D.*G;
    end
    S(ii) = abs(P);
end


% ===== Plotting =========================================================
if(useplot)
    figure; semilogx(f,db(S));
    ylabel('Amplitude (dB)');
    xlabel('Frequency (Hz)');
end
