function conf = SFS_config()
%SFS_CONFIG Configuration file for the SoundFieldSynthesis functions
%
%   Usage: conf = SFS_config
%
%   Output parameters:
%       conf    - struct containing all configuration variables
%
%   SFS_CONFIG() creates the struct conf containing the default
%   configuration values. If you want to create other entries, please set
%   them in your script (e.g. conf.fs = 48000) and pass the conf struct to
%   the desired function as last input (e.g. tapering_window(x0,conf)).
%
%   So edit this function only, if the default values have changed!
%
%   see also: SFS_start

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 0;
nargmax = 0;
error(nargchk(nargmin,nargmax,nargin));


%% ===== Configuration default values ===================================

% ===== Table of Content ========================
%
% - Misc
% - Audio
% - Simulations
% - Secondary Sources
% - Binaural Reproduction
%   * Headphone compensation
%   * HRIR/BRIR
%   * Auralization
% - WFS
%   * Pre-equalization
%   * Tapering
%   * Virtual Sources
% - SDM
% - HOA
% - Plotting
%   * Gnuplot
%


% ===== Misc ====================================
conf.tmpdir = '/tmp/sfs';
% Debugging level. We are supporting 3 levels:
%   0 - normal mode
%   1 - debug modus, showing interim results and plots
conf.debug = 0;


% ===== Audio ===================================
% Samplingrate
conf.fs = 44100; % Hz
% Speed of sound
conf.c = 343; % m/s
% use fractional delays for delay lines
conf.usefracdelay = 0;
conf.fracdelay_method = 'resample';


% ===== Simulations =============================
% xyresolution for wavefield simulations
conf.xysamples = 300; % samples
% Phase of omega of wavefield (change this value to create monochromatic wave
% fields with different phases for a movie)
conf.phase = 0; % rad
% Time frame to simulate for wave field in time domain (change this value to
% create impulse response movies of the wave field
conf.frame = 1000;


% ===== Secondary Sources =======================
% Interspacing (distance) between the secondary sources
conf.dx0 = 0.15; % m
% Array position
conf.X0 = [0 0 0]; % m
% Array geometry
% Possible values are: 'linear', 'box', 'circle', 'U', 'custom'
conf.array = 'linear';
% Vector containing custom secondary source positions and directions.
% conf.x0 = [x0; y0; z0; phi];
conf.x0 = []; % m; m; m; rad


% ===== Binaural reproduction ===================
%
% === Headphone compensation ===
% Headphone compensation
conf.usehcomp = true; % boolean
% Headphone compensation file for left and right ear.
conf.hcomplfile = ...
    '~/data/ir_databases/headphone_compensations/TU_FABIAN_AKGK601_hcompl.wav';
conf.hcomprfile = ...
    '~/data/ir_databases/headphone_compensations/TU_FABIAN_AKGK601_hcompr.wav';
%
% === HRIR/BRIR ===
% Target length of BRIR impulse responses (2^14 may be enough for your
% purposes, but for a large distance between source and listener, this will
% be not enough to contain the desired time delay. But don't worry, SFS
% checks for you if conf.N is large enough)
conf.N = 2^15; % samples
% To use a dynamic binaural simulation together with the SoundScape Renderer
% (SSR) and a headtracker, brs sets can be created. If these sets should be
% used in BRS mode of the SSR, the angles have to be:
% conf.brsangles = 0:1:359;
% If the brs set should be used as IRs for the SSR, the angles have to be:
% conf.brsangles = 360:-1:1;
conf.brsangles = 0:1:359; % degree
%
% === Auralisation ===
% These files are used for the auralization of impulse responses by the
% auralize_ir() function.
% NOTE: you have to provide them by yourself!
conf.speechfile = '';
conf.cellofile = '';
conf.castanetsfile = '';
conf.noisefile = '';
conf.pinknoisefile = '';


% ===== WFS =====================================
% The amplitude will be correct at the point xref for 2.5D
% synthesis.
% Thi point is also used to scale the wave field to 1 at this point.
conf.xref = [0 2 0]; % m, m, m
%
% ===== Pre-Equalization =====
% WFS can be implemented very efficiently using a delay-line with different
% amplitudes and convolving the whole signal once with the so called
% pre-equalization filter [References]. If we have aliasing in our system we
% only want to use the pre-equalization filter until the aliasing frequency,
% because of the energy the aliasing is adding to the spectrum above this
% frequency (which means the frequency response over the aliasing frequency is
% allready "correct") [Reference]
% Use WFS preequalization-filter
conf.usehpre = false; % boolean
% Lower frequency limit of preequalization filter (~ frequency when
% subwoofer is active)
conf.hpreflow = 50; % Hz
% Upper frequency limit of preequalization filter (~ aliasing frequency of
% system)
conf.hprefhigh = 1200; % Hz
%
% ===== Tapering =====
% The truncation of the loudspeaker array leads to diffraction of the
% synthesized wave field. It has been shown that the truncation can be discribed
% by cylindrical waves originating from the edges of the array
% [Young,Sommerfeld,Rubinovitch]. Therefore a good method to reduce artifacts
% due to the diffraction edge waves is to fade out the amplitude of the driving
% function at the edges of the array. This method is called tapering and
% implemented using a Hanning window.
% Use tapering window
conf.usetapwin = true; % boolean
% Size of the tapering window
conf.tapwinlen = 0.3; % percent of array length 0..1
%
% === Virtual Sources ===
% Pre-delay for causality for focused sources
% Note: for a non focused point source this will be set automaticly to 0
% from brs_wfs_25d!
conf.t0 = -3*1024/conf.fs; % s


% ===== SDM =====================================
% Use the evanescent part of the driving function for SDM
conf.withev = true; % boolean


% ===== HOA =====================================
% TODO


% ===== Plotting ================================
% Plot the results (wave fields etc.) directly
conf.useplot = false; % boolean
% Plot mode (uses the GraphDefaults function). Avaiable modes are:
%   'monitor'   - displays the plot on the monitor
%   'paper'     - eps output in conf.plot.outfile
%   'png'       - png output in conf.plot.outfile
conf.plot.mode = 'monitor';
% Plot amplitudes in dB (e.g. wavefield plots)
conf.plot.usedb = false; % boolean
% caxis settings (leave blank, if you would use the default values of the given
% plot function)
conf.plot.caxis = '';
% Default colormap to use (note this works only with Matlab at the moment)
conf.plot.colormap = 'gray'; % (default: 'gray')
% Plot loudspeakers in the wave field plots
conf.plot.loudspeakers = true; % boolean
% Use real loudspeakers symbols (otherwise crosses are used)
conf.plot.realloudspeakers = true; % boolean
% Size of the loudspeaker
conf.plot.lssize = conf.dx0; % m
% Size of the plot
conf.plot.size = [16,11.55]; % cm
% Additional plot command
conf.plot.cmd = '';
% output of plot (file or screen)
conf.plot.usefile = 0;
conf.plot.file = 'sfs';
%
% === Gnuplot ===
% Use gnuplot
conf.plot.usegnuplot = false; % boolean
