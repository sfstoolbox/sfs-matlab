%AES127_ECHO_DIRECTION Calculates the direction of the pre-echos for the
%positions used in AES127 paper
%
% L = 30m, 10m; (xs,ys) = (0,1)
% (X,Y) = (0,3), (2,4), (5,4), (10,4), (15,4); phi = 0 for all positions
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
% Listener direction offset (defines the 0° direction of the listener,
% default: 0° == negative y-direction)
% This value is the reference direction for every angle value given to the
% SFS functions (so if you change this value to 90, your other angles have 
% to change -90).
conf.listoffset = 0;
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

% === 10m Array ===
L = 10;
X = [0 2 5];
Y = [3 4 4];
for i = 1:3
    echo_direction(X(i),Y(i),0,xs,ys,L,conf);
end

% === 30m Array ===
L = 30;
X = [0 2 5 10 15];
Y = [3 4 4 4 4];
for i = 1:5
    echo_direction(X(i),Y(i),0,xs,ys,L,conf);
end
