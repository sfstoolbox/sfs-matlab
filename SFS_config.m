function conf = SFS_config()
%SFS_CONFIG Configuration file for the SoundFieldSynthesis functions
%   Usage: conf = SFS_config
%
%   Output parameters:
%       conf    - struct containing all configuration variables
%
%   SFS_CONFIG() creates the struct conf containing the default
%   configuration values. If you want to create other entries, please set
%   them in your script (e.g. conf.fs = 48000) and pass the conf struct to
%   the desired function as last input (e.g. tapwin(L,conf)).
%
%   So edit this function only, if the default values have changed!
%
%   see also:
%

%   FIXME:
%       * can the hpre.fhigh (aliasing frequency) be calculated, or should
%         it be set by hand? (Bug #10)
%         UPDATE: one should be able to calculate this using diffraction theory

% AUTHOR: Hagen Wierstorf

%% ===== Checking of input  parameters ==================================

if nargchk(0,0,nargin)
    error('Wrong number of args. Usage: conf = SFS_config');
end


%% ===== Configuration default values ===================================

% ===== Path settings ===========================
conf.sfspath = '/home/hagen/svn/sfs';
conf.tmpdir = '/tmp/sfs';


% ===== General audio ===========================
% Samplingrate (Hz)
conf.fs = 44100;
% Speed of sound (m/s)
conf.c = 343;


% ===== Binaural reproduction ===================
% Headphone compensation (true or false)
conf.usehcomp = true;
% Headphone compensation file for left and right ear.
conf.hcomplfile = ...
    ['~/data/measurements/headphone_compensation_filters/'...
     'FABIAN_AKG_K601/'...
     'FABIAN_AKG_K601_(left)_4096_min_phase_marg100_inverse.wav'];
conf.hcomprfile = ...
    ['~/data/measurements/headphone_compensation_filters/'...
     'FABIAN_AKG_K601/'...
     'FABIAN_AKG_K601_(right)_4096_min_phase_marg100_inverse.wav'];


% ===== WFS =====================================
% Loudspeaker distance (m)
conf.LSdist = 0.15;
% Array position (m)
conf.X0 = 0;
conf.Y0 = 0;
% Array geometry
% Possible values are: 'linear', 'box', 'circle', 'U', 'custom'
conf.array = 'linear';
% Vector containing custom loudspeaker positions
conf.x0 = [];
conf.y0 = [];
conf.phi = [];
% Listener direction offset (defines the 0° direction of the listener,
% default: 0° == negative y-direction)
% This value is the reference direction for every angle value given to the
% SFS functions (so if you change this value to 90, your other angles have
% to change -90).
conf.listoffset = 0;
% y position for amplitude calculation
% NOTE: the amplitude will be correct at the line parallel to the x-axis given
% at the y position (this is a 2.5D property of WFS)
conf.yref = 2;
%
% ===== Pre-Equalization =====
%
% WFS can be implemented very efficiently using a delay-line with different
% amplitudes and convolving the whole signal once with the so called
% pre-equalization filter [References]. If we have aliasing in our system we
% only want to use the pre-equalization filter until the aliasing frequency,
% because of the energy the aliasing is adding to the spectrum above this
% frequency (which means the frequency response over the aliasing frequency is
% allready "correct") [Reference]
%
% WFS preequalization-filter (true or false)
conf.usehpre = false;
% Lower frequency limit of preequalization filter (= frequency when
% subwoofer is active) (Hz)
conf.hpreflow = 50;
% Upper frequency limit of preequalization filter (= aliasing frequency of
% system) (Hz)
conf.hprefhigh = 1200;
%
% ===== Tapering =====
%
% The truncation of the loudspeaker array leads to diffraction of the
% synthesized wave field. It has been shown that the truncation can be discribed
% by cylindrical waves originating from the edges of the array
% [Young,Sommerfeld,Rubinovitch]. Therefore a good method to reduce artifacts
% due to the diffraction edge waves is to fade out the amplitude of the driving
% function at the edges of the array. This method is called tapering and
% implemented using a Hanning window.
%
% Use tapering window (true or false)
conf.usetapwin = true;
% Size of the tapering window (% of array length => 0..1)
conf.tapwinlen = 0.3;


% ===== SDM =====================================
% Use the evanescent part of the driving function for SDM (true or false)
conf.withev = true;


% ===== HOA =====================================


% === Virtual Source ===
% Pre-delay for causality for focused sources (s)
% Note: for a non focused point source this will be set automaticly to 0
% from wfs_brs!
conf.t0 = -3*1024/conf.fs;


% === HRIR/BRIR ===
% IR dataset
conf.irsfile = '~/data/measurements/HRIRs/FABIAN_postprocessed_anechoic.mat';
% Target length of BRIR impulse responses (2^14 may be enough for your
% purposes, but for a large distance between source and listener, this will
% be not enough to contain the desired time delay. But don't worry, SFS
% checks for you if conf.N is large enough)
conf.N = 2^15;
% Angles for the BRS set for the SSR (first two columns of the BRIR are 0° of
% the listener direction, next two columns are -1°, etc. So the brsangles has

% to spin around the other direction, because the source moves in opposite
% direction to the listener direction; see wfs_brs_set)
% NOTE: This has changed recently. Our old approach uses a angle notation that
% has had a counterwise rotation for positive values. The new one uses a
% counterclockwise notation (as in Mathematics) for positive values.
% The old brsangles notation was therefore: conf.brsangles = 360:-1:1;
conf.brsangles = 0:1:359;
% Auralisation files
conf.speechfile = '~/data/signals/goesa_sentence.wav';
conf.cellofile = '~/data/signals/cello.wav';
conf.castanetsfile = '~/data/signals/castanets.wav';
conf.noisefile = '~/data/signals/noise.wav';
conf.pinknoisefile = '~/data/signals/pinknoise.wav';


% === Wave Field Simulation ===
% xyresolution for wavefield simulations
conf.xysamples = 300;
% Phase of omega (rad)
conf.phase = 0;
% Time frame to simulate for wave field in time domain
conf.frame = 1000;


% === Plotting ===
% Plot the results (wave fields etc.) directly
conf.useplot = false;
% Use gnuplot
conf.plot.usegnuplot = false;
% Plot mode (uses the GraphDefaults function). Avaiable modes are:
%   'monitor'   - displays the plot on the monitor
%   'paper'     - eps output in conf.plot.outfile
%   'png'       - png output in conf.plot.outfile
conf.plot.mode = 'monitor';
% Outfile for the 'paper' and 'png' plot modes
conf.plot.outfile = 'sfs';
% Plot amplitudes in dB (e.g. wavefield plots)
conf.plot.usedb = false;
% caxis settings (leave blank, if you would use the default values of the given
% plot function)
conf.plot.caxis = '';
% Plot loudspeakers in the wave field plots
conf.plot.loudspeakers = true;
% Use real loudspeakers symbols (otherwise crosses are used)
conf.plot.realloudspeakers = false;
% Size of the loudspeaker
conf.plot.lssize = 0.1;
