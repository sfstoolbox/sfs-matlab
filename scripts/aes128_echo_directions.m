%AES128_ECHO_DIRECTION Calculates the direction of the pre-echos
%
%

% AUTHOR: Hagen Wierstorf


%% ===== Configuration ==================================================

% === General audio ===
% Samplingrate (Hz)
conf.fs = 44100;
% Speed of sound (m/s)
conf.c = 343;


% === WFS === 
% Loudspeaker distance (m)
conf.dx0 = 0.15;
% Array position (m)
conf.X0 = 0;                    
conf.Y0 = 0;
% WFS preequalization-filter (true or false)
conf.usehpre = false;
% Lower frequency limit of preequalization filter (= frequency when 
% subwoofer is active) (Hz)
conf.hpreflow = 50;
% Upper frequency limit of preequalization filter (= aliasing frequency of 
% system) (Hz)
conf.hprefhigh = 1200;
% Use tapering window
conf.usetapwin = true;
% Size of the tapering window (% of array length => 0..1)
conf.tapwinlen = 0.3;


% === Virtual Source ===
% Pre-delay for causality for focused sources (s)
% Note: for a non focused point source this will be set automaticly to 0
% from wfs_brs!
conf.t0 = -3*1024/conf.fs;


% === Plotting ===
% Plot the results etc.
conf.useplot = false;
% Use gnuplot
conf.usegnuplot = true;


%% ===== Computation ====================================================
% Focused source position
xs = 0;
ys = 1;

% === 4m Array ===
L = 4;
R = 1;
% Listener direction (0°,-30°,-60°)
for alpha = [0 -30 -60]
    X = R * sin(-alpha/180*pi) + xs;
    Y = sqrt(R^2-(X-xs)^2) + ys;
    echo_direction(X,Y,alpha,xs,ys,L,conf);
end

% === 10m Array ===
L = 10;
R = 4;
% Listener direction (0°,-30°,-60°)
for alpha = [0 -30 -60]
    X = R * sin(-alpha/180*pi) + xs;
    Y = sqrt(R^2-(X-xs)^2) + ys;
    echo_direction(X,Y,alpha,xs,ys,L,conf);
end

