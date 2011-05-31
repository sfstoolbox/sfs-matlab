% Create BRS sets for a linear array

% The data of the corresponding experiment is located here:
% https://www.qu.t-labs.tu-berlin.de/audio/experiments/2010-06_AES129_focused_sources
%
% The HRIR dataset measured with FABIAN is located here:
% https://www.qu.t-labs.tu-berlin.de/audio/measurements/HRIRs/FABIAN_anechoic
%
% The headphone compensation can be found here:
% https://www.qu.t-labs.tu-berlin.de/audio/measurements/headphone_compensation_filters/FABIAN_AKG_K601
% 

% AUTHOR: Hagen Wierstorf


%% ===== Variables ======================================================

clear all;

% Geometry of the test
%
%    x-axis                      [X0 Y0]
%       <------^--^--^--^--^--^--^--^--^--^--^--^--^--^--^-------
%                                   |
%                                   x [xs ys] (Focused source)
%                                 | |
%                               |   |       
%                             |     |       
%                           | alpha |  
%                         | R       |
%                       |           |
%                     O             |
%                  [X,Y]            |
%                                   v y-axis
%
% Radius of listener positions
R = 1;
% Angle of listener position (Note: if no angle is given, only the
% reference BRS set using a single source is created, otherwise also BRS
% sets for focused sources located at the given direction are created
alpha = [0 -30 -60];
% Aliasing frequency for the given position
fal = [5000 4100 2200];     % L = 4, R = 1;
%fal = [3700 2100 1500];     % L = 10; R = 4;
%fal = [1500 1500 1500];
% Array length
%L = [0.5 1 2 10];
L = [0.5 1 2 4];
% Focused source position
xs = 0;
ys = 1;
% Directory to store BRS sets
outdir = ['~/d/data/experiments/2010-06_AES129_focused_sources/'...
          'brirs/'];


%% ===== Configuration ==================================================

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
    ['~/d/data/measurements/headphone_compensation_filters/'...
     'FABIAN_AKG_K601/'...
     'FABIAN_AKG_K601_(left)_4096_min_phase_marg100_inverse.wav'];
conf.hcomprfile = ...
    ['~/d/data/measurements/headphone_compensation_filters/'...
     'FABIAN_AKG_K601/'...
     'FABIAN_AKG_K601_(right)_4096_min_phase_marg100_inverse.wav'];
 

% === WFS === 
% Loudspeaker distance (m)
conf.dx0 = 0.15;
% Array position (m)
conf.X0 = 0;                    
conf.Y0 = 0;
% WFS preequalization-filter (true or false)
conf.usehpre = true;
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


% === HRIR ===
% HRIR dataset
conf.irsfile = '~/d/data/measurements/HRIRs/FABIAN_postprocessed_anechoic.mat';


% === BRIR ===
% Target length of BRIR impulse responses
conf.N = 2^14;
% Angles for the BRS set for the SSR (first two columns of the BRIR are 0° of 
% the listener direction, next two columns are -1°, etc. So the brsangles has 
% to spin around the other direction, because the source moves in opposite 
% direction to the listener direction; see wfs_brs_set)
conf.brsangles = 0:1:359;
% Auralisation files
conf.speechfile = '~/d/data/signals/speech.wav';
conf.cellofile = '~/d/data/signals/cello.wav';
conf.castanetsfile = '~/d/data/signals/castanets.wav';


% === Plotting ===
% Plot the results etc.
conf.useplot = false;
% Use gnuplot
conf.usegnuplot = false;



%% ===== Computation ====================================================

% Read HRIRs
hrirs = read_irs(conf);

% Compute BRIR for reference single source
brir = ref_brs_set(0,R+1,0,xs,ys,hrirs,conf);

% Scale BRIR output (|BRIR|<1)
%brir = 0.95*brir/max(abs(brir(:)));
brir = brir ./ max([rms(brir(:,1)) rms(brir(:,2))]) .* 0.001;
% Write BRIR to wave file
outfile = sprintf('%sbrir_ref_xs%d_ys%d_R%d_phi0.wav',outdir,xs,ys,R);
wavwrite(brir,conf.fs,16,outfile);

% Compute BRIRs
%for l = 1:length(L)
%    
%    % Different listener directions
%    for a = 1:length(alpha)
%        
%        % Calculate listener positions (X,Y)
%        X = R * sin(-alpha(a)/180*pi) + xs;
%        Y = sqrt(R^2-(X-xs)^2) + ys;
%        
%        % Set fhigh for the WFS pre-equalization filter (aliasing
%        % frequency)
%        conf.hprefhigh = fal(a);
%        
%        % Print (and set) outfile, so we know were the script is
%        if(conf.usehpre)
%            outfile = ...
%                sprintf('%sbrir_%dm_xs%d_ys%d_R%d_phi%d_fal%d.wav',...
%                        outdir,L(l),xs,ys,R,alpha(a),fal(a));
%        else
%            outfile = sprintf('%sbrir_%dm_xs%d_ys%d_R%d_phi%d.wav',...
%                              outdir,L(l),xs,ys,R,alpha(a));
%        end
%        
%        % Calculate BRIR
%        brir = wfs_brs_set(X,Y,alpha(a),xs,ys,L(l),hrirs,conf);
%        
%        % Scale BRIR output (|BRIR|<1) using the first channel (0°)
%        brir = brir ./ max([rms(brir(:,1)) rms(brir(:,2))]) .* 0.001;
%        %brir = 0.95*brir/max(abs(brir(:)));
%        % Write BRIR to wave file
%        wavwrite(brir,conf.fs,16,outfile);
%    end
%end
