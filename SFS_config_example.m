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
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
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
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 0;
nargmax = 0;
narginchk(nargmin,nargmax);


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


%% ===== Misc ============================================================
conf.tmpdir = '/tmp/sfs'; % string
% Debugging level. We are supporting 3 levels:
%   0 - normal mode
%   1 - debug modus, showing interim results and plots
conf.debug = 0; % 0 or 1
% Show a progress bar in the loops (for example sound_field_mono). This can be
% useful if you are using secondary sources with >1000 loudspeakers
conf.showprogress = false; % boolean


%% ===== Audio ===========================================================
% Audio and signal processing settings
%
% Samplingrate
conf.fs = 44100; % / Hz
% Speed of sound
conf.c = 343; % / m/s
% use fractional delays for delay lines
conf.usefracdelay = false; % boolean
conf.fracdelay_method = 'resample'; % string
% Bandpass filter applied in sound_field_imp()
conf.usebandpass = true; % boolean
conf.bandpassflow = 10; % / Hz
conf.bandpassfhigh = 20000; % / Hz


%% ===== SFS =============================================================
% Common sound field synthesis settings
%
% === Dimensionality ===
% Dimensionality of the secondary sources and the sound field synthesis driving
% functions:
% '2D'    - line sources as secondary sources, arranged in a circle, line, ...
% '2.5D'  - point sources as secondary sources, arranged in a circle, line, ...
% '3D'    - point sources as secondary sources, arranged in a sphere, plane, ...
conf.dimension = '2.5D'; % string
%
% === Driving functions ===
% Implementation of driving functions. For the default ones use 'default'. These
% functions are described in the PDF documentation, in the doc folder of the
% SFS-Toolbox. For possible other flags have a look into the driving functions.
% Most users can safely use the 'default' flag here.
conf.driving_functions = 'default'; % string
%
% === Impulse responses ===
% Length of impulse responses used in the time domain driving functions
% and for the creation of the binaural simulations.
% Don't worry, SFS checks for you if conf.N is large enough)
conf.N = 2048; % samples
%
% === 2.5D ===
% The amplitude will be correct at the point xref for 2.5D
% synthesis.
% This point is also used to scale the sound field to 1 at this point.
conf.xref = [0 0 0]; % / m
%
% === Tapering ===
% The truncation of the loudspeaker array leads to diffraction of the
% synthesized sound field. It has been shown that the truncation can be discribed
% by cylindrical waves originating from the edges of the array
% [Young,Sommerfeld,Rubinovitch]. Therefore a good method to reduce artifacts
% due to the diffraction edge waves is to fade out the amplitude of the driving
% function at the edges of the array. This method is called tapering and
% implemented using a Hanning window.
% Use tapering window
conf.usetapwin = true; % boolean
% Size of the tapering window
conf.tapwinlen = 0.3; % / percent of array length, 0..1


%% ===== Wave Field Simulations ==========================================
% Simulations of monochromatic or time domain sound field
%
% xyz-resolution for sound field simulations, this value is applied along every
% desired dimension, expect if only one point is desired
conf.resolution = 300; % / samples
% Phase of omega of sound field (change this value to create monochromatic sound
% fields with different phases, for example this can be useful to create a movie)
conf.phase = 0; % / rad
% normalize the simulated sound field?
conf.usenormalisation = true; % boolean


% ===== Secondary Sources =======================
% Number of secondary sources
conf.secondary_sources.number = 64; % integer
% Diameter/Length of secondary source array
conf.secondary_sources.size = 3; % / m
% Center of array, X0
conf.secondary_sources.center = [0 0 0]; % / m
% Array geometry
% Possible values are: 'line', 'box', 'circle', 'sphere'
conf.secondary_sources.geometry = 'circle'; % string
% Vector containing custom secondary source positions and directions.
% conf.secondary_sources.x0 = [x0; y0; z0; nx0; ny0; nz0];
conf.secondary_sources.x0 = []; % / m
% Grid for the spherical array. Note, that you have to download and install the
% spherical grids from an additiona source. For available grids see:
% http://github.com/sfstoolbox/data/tree/master/spherical_grids
% An exception are Gauss grids, which are available via 'gauss' and will be
% calculated on the fly allowing very high number of secondary sources.
conf.secondary_sources.grid = 'equally_spaced_points'; % string


%% ===== WFS =============================================================
% Settings for wave field synthesis
%
% === Pre-Equalization ===
% WFS can be implemented very efficiently using a delay-line with different
% amplitudes and convolving the whole signal once with the so called
% pre-equalization filter [References]. If we have aliasing in our system we
% only want to use the pre-equalization filter until the aliasing frequency,
% because of the energy the aliasing is adding to the spectrum above this
% frequency (which means the frequency response over the aliasing frequency is
% allready "correct") [Reference]
% Use WFS preequalization-filter
conf.wfs.usehpre = true; % boolean
% FIR or IIR pre-equalization filter
% NOTE: only FIR is working under octave at the moment
conf.wfs.hpretype = 'FIR'; % 'FIR' or 'IIR'
% Lower frequency limit of preequalization filter (~ frequency when
% subwoofer is active)
conf.wfs.hpreflow = 50; % / Hz
% Upper frequency limit of preequalization filter (~ aliasing frequency of
% system)
conf.wfs.hprefhigh = 1200; % / Hz
% IIR bandwidth for the Lagrange interpolation region
conf.wfs.hpreBandwidth_in_Oct = 2; % / octaves
% desired IIR filter order
conf.wfs.hpreIIRorder = 4; % integer


%% ===== SDM =============================================================
% Use the evanescent part of the driving function for SDM
conf.sdm.withev = true; % boolean


%% ===== HOA =============================================================
% normally the order of NFC-HOA is set by the nfchoa_order() function which
% returns the highest order for which no aliasing occurs. If you wish to use
% another order you can set it manually here, otherwise leave it blank
conf.nfchoa.order = []; % integer


%% ===== Binaural reproduction ===========================================
% Settings regarding all the stuff with impulse responses from the SFS_ir and
% SFS_binaural_synthesis folders
%
% Directory containing HRTF data bases, you want to use. Note, that also all
% subdirectories will be added to the path. This is not done automatically, but
% by calling addirspath;
% If you have more than one path, seperate them by :
conf.ir.path = '~/git/sfs/data/HRTFs:~/svn/ir_databases:~/svn/measurements'; % string
%
% If we load an HRTF data set we are most likely interested to modify its
% existing length, to enable a delaying of the impulse responses without
% problems. If these value is set to "false", zeros are padded at the
% beginning of all HRTFs corresponding to the maximum distance of the whole
% set. In addition the overall length of the impulse responses is set to
% conf.N. This is applied directly if you load a HRTF set with read_irs().
% Only set this to true if you really know what you are doing.
conf.ir.useoriglength = false; % boolean
%
% Use interpolation to get the desired HRTF for binaural simulation. If this is
% disabled the HRTF returned by a nearest neighbour search is used instead.
% Depending on the geometry of the measured HRTF data set, the interpolation
% will be done between the two or three nearest HRTFs.
conf.ir.useinterpolation = true; % boolean
%
% === Headphone compensation ===
% Headphone compensation
conf.ir.usehcomp = true; % boolean
% Headphone compensation file for left and right ear.
conf.ir.hcompfile = 'data/headphone_compensation/QU_KEMAR_AKGK601_hcomp.wav'; % string
%
% === Auralisation ===
% These files are used for the auralization of impulse responses by the
% auralize_ir() function.
% NOTE: you have to provide them by yourself!
conf.ir.speechfile = ''; % string
conf.ir.cellofile = ''; % string
conf.ir.castanetsfile = ''; % string
conf.ir.noisefile = ''; % string
conf.ir.pinknoisefile = ''; % string
%
% === SoundScape Renderer ===
% To use a dynamic binaural simulation together with the SoundScape Renderer
% (SSR) and a headtracker, brs sets can be created. If these sets should be
% used in BRS mode of the SSR, the angles have to be:
% conf.ir.brsangles = 0:1:359;
% If the brs set should be used as IRs for the SSR, the angles have to be:
% conf.ir.brsangles = 360:-1:1;
conf.ir.brsangles = 0:1:359; % / degree


%% ===== Plotting ========================================================
% Plot the results (sound fields etc.) directly
conf.plot.useplot = false; % boolean
% Plot mode (uses the GraphDefaults function). Avaiable modes are:
%   'monitor'   - displays the plot on the monitor
%   'paper'     - eps output in conf.plot.outfile
%   'png'       - png output in conf.plot.outfile
conf.plot.mode = 'monitor'; % string
% Plot amplitudes in dB (e.g. sound field plots)
conf.plot.usedb = false; % boolean
% caxis settings (leave blank, if you would use the default values of the given
% plot function)
conf.plot.caxis = []; % [min max]
% Default colormap to use
% The Toolbox comes with two own color maps, if you set 'default' or 'moreland'
% you will get a blue/red-colormap after
% http://www.sandia.gov/~kmorel/documents/ColorMaps/
% If you set 'gray' or 'grey' you will get a colormap ranging from white to
% black. In addition you can add every other map you can specify in
% Matlab/Octave. For example to get the Matlab default colormap ser 'jet'.
conf.plot.colormap = 'default'; % string
% Plot loudspeakers in the sound field plots
conf.plot.loudspeakers = true; % boolean
% Use real loudspeakers symbols (otherwise crosses are used)
conf.plot.realloudspeakers = false; % boolean
% Size of the loudspeaker
conf.plot.lssize = 0.16; % m
% Size of the plot
conf.plot.size_unit = 'px'; % 'px','cm','inches'
conf.plot.size = [540 404]; % [xsize ysize]
% Resolution of plot in dpi
conf.plot.resolution = 150; % integer
% Additional plot command
conf.plot.cmd = ''; % string
% output of plot (file or screen)
conf.plot.usefile = false; % boolean
% File name, if this is provided with as *.png or *.eps file, the figure is
% plotted to the regarding file
conf.plot.file = ''; % string
%
% === Gnuplot ===
% Use gnuplot
conf.plot.usegnuplot = false; % boolean
