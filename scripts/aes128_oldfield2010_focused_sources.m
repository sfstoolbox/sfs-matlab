% Analysing the results of the paper:
% R. Oldfield, U. Drumm, J. Hirst, "The Perception of Focused Sources in
% Wave Field Synthesis as a Function of Listener Angle", AES128th,
% May 2010, London.
%
% In this paper the localization for different focused source position in a
% rectengular array was investigated. 
%
% Geometry:
%
%     |
%     |             3 m, 22 loudspeakers
%     # # # # # # # # # # # # # # # # # # # # # #
%     #                                         #
%     #                                         #
%     #                                         #
%     #                                       --#
%     #                                    ---  #
%     #                                 ---     #
%     #                              ---        #
%     #                           ---           #
%     #                   x  x ---              #
%     #               x     ---  x (xs,ys)      #
%     #---               ---                    #
%     #   ---     x   ---           x           #
%     #      ---   --- \                        # 4 m, 30 loudspeaker
%     #----------x--phi-|--O-(X,Y)----x---------#-----> y
%     #      ---   --- /                        #
%     #   ---     x   ---           x           #
%     #---               ---                    #
%     #               x     --- x               #
%     #                   x  x ---              #
%     #                           ---           #
%     #                              ---        #
%     #                                 ---     #
%     #                                    ---  #
%     #                                       --#
%     #                                         #
%     #                                         #
%     #                                         #
%     #                                         #
%     # # # # # # # # # # # # # # # # # # # # # #
%     |
%     |
%     \/
%     x
%
% Parameter of his array:
% 3 m x 4 m, 104 loudspeaker => 22, 30
% => loudspeaker distance = 0.135 m
%
% View angle:
% phi = 70° constant for all focused source positions => different number
% of loudspeakers active per position of focused source
%

% AUTHOR: Hagen Wierstorf


%% ===== Variables ======================================================

clear all;

% Radius of listener positions: 0.45, 0.6, 0.75, 0.9 m
%R = [0.45,0.6,0.75,0.9];
R = 0.9;
% View angle
phi = 70/180*pi;
% Focused source angles
psi_fs = [0:-15:-180,-187.5:-15:-352.5]./180*pi;
% Used coordinates
x_coord = -2:0.1:2;
y_coord = -1.5:0.1:1.5;
% Listener position
X = 0;
Y = 0;

% Monochromatic frequency
f = 2000;

% Content to produce (castanets, speech, cello)
content = 'castanets';
% Directory to store the wave files
%outdir = ['D:/data/experiments/2010-01_AES128_focused_sources/'...
%          'stimuli/enhanced/'];
outdir = 'D:\data\experiments\2010-01_AES128_focused_sources\stimuli\enhanced\';


%% ===== Configuration ==================================================

% Get default values
conf = SFS_config;

% Explicit configuration to ensure reproducibility
% === General audio ===
% Samplingrate (Hz)
conf.fs = 44100;
% Speed of sound (m/s)
conf.c = 343;


% === Headphone ===
% Headphone compensation (true or false)
conf.usehcomp = true;
% Headphone compensation file for left and right ear.
conf.hcomplfile = ...
    ['D:/data/measurements/headphone_compensation_filters/'...
     'FABIAN_AKG_K601/'...
     'FABIAN_AKG_K601_(left)_4096_min_phase_marg100_inverse.wav'];
conf.hcomprfile = ...
    ['D:/data/measurements/headphone_compensation_filters/'...
     'FABIAN_AKG_K601/'...
     'FABIAN_AKG_K601_(right)_4096_min_phase_marg100_inverse.wav'];


% === WFS === 
% Loudspeaker distance (m)
conf.LSdist = 0.135;
% Array position (m)
conf.X0 = 0;                    
conf.Y0 = -1.5;
% Listener direction offset (defines the 0° direction of the listener,
% default: 0° == negative y-direction)
conf.listoffset = 0;
conf.yref = 1.3;
% WFS preequalization-filter (true or false)
conf.usehpre = true;
% Lower frequency limit of preequalization filter (= frequency when 
% subwoofer is active) (Hz)
conf.hpreflow = 25;
% Upper frequency limit of preequalization filter (= aliasing frequency of 
% system) (Hz)
conf.hprefhigh = 2500;
% Use tapering window
conf.usetapwin = true;
% Size of the tapering window (% of array length => 0..1)
conf.tapwinlen = 0.3;
% Use the evanescent part of the driving function (SDM method, see ...)
conf.withev = 0;


% === Virtual Source ===
% Pre-delay for causality for focused sources (s)
% Note: for a non focused point source this will be set automaticly to 0
% from wfs_brs!
conf.t0 = -3*1024/conf.fs;


% === HRIR ===
% HRIR dataset
conf.irsfile = '~/d/data/measurements/HRIRs/FABIAN_postprocessed_anechoic.mat';


% === BRIR ===
% Target length of BRIR impulse responses
conf.N = 2^14;
% Angles for the BRS set for the SSR (first two columns of the BRIR are 0° 
% of the listener direction, next two columns are 1°, etc. So the 
% brsangles has to spin around the other direction, because the source 
% moves in opposite direction to the listener direction; see wfs_brs_set)
conf.brsangles = 360:-1:1;
% Auralisation files
conf.speechfile = 'D:/data/signals/speech.wav';
conf.cellofile = 'D:/data/signals/cello.wav';
conf.castanetsfile = 'D:/data/signals/castanets.wav';


% === Plotting ===
% Plot the results etc.
conf.useplot = false;
% Use gnuplot
conf.usegnuplot = false;
% Plot amplitudes in dB (e.g. wavefield plots)
conf.usedb = 0;


%% ===== Computation ====================================================

echo_time_listener = zeros(length(psi_fs),length(R));
cond = 1;
for jj = 1:length(R)
    for ii = 1:length(psi_fs)
    
        % Calculate focused source position
        ys = -R(jj)*cos(psi_fs(ii));
        xs =  R(jj)*sin(psi_fs(ii));
    
        % Calculate array length (which means active loudspeakers)
        % L/2 = |y_min-ys| * tan(phi/2)
        L = 2 * abs(conf.Y0-ys) * tan(phi/2);
        
        % Set the x center of the loudspeaker array to the focused source
        % position
        conf.X0 = xs;
        
        % Calculate and plot the wavefield for the given focused source
        % position
        [x,y,P] = SDM_25D_wavefield(4,3,xs,ys,L,f,conf);
        plot_wavefield(x,y,P,conf); hold on;
        
        % Draw a circle as head
        d = 0.1;
        [x,y,z] = cylinder(d/2,200);
        plot(x(1,:),y(1,:),'-w','LineWidth',12);
        title_str = sprintf('Condition %2.0f: %.0f°; [xs,ys] = %.2f, %.2f; L = %.2f',...
            cond,psi_fs(ii)/pi*180,xs,ys,L);
        title(title_str);
        
        %plot_echo_times(x_coord,y_coord,xs,ys,L,conf);
        
        % Calculate the echo times for the given focused source position
        echo_time_listener(ii,jj) = echo_time(X,Y,xs,ys,L,conf);
        
        
        cond = cond+1;
        
    end
end


figure;
plot(psi_fs/pi*180,echo_time_listener(:,1)*1000,'-b'); hold on;
%plot(psi_fs/pi*180,echo_time_listener(:,2)*1000,'-g'); hold on;
%plot(psi_fs/pi*180,echo_time_listener(:,3)*1000,'-r'); hold on;
%plot(psi_fs/pi*180,echo_time_listener(:,4)*1000,'-m');
xlabel('Direction (°)');
ylabel('t (ms)');
axis([-360,0,-6,0]);
title('Arriving time of first pre-echo');
