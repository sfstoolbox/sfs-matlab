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
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
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
% - Delayline
% - Sound Field Synthesis (SFS)
%   * Dimensionality
%   * Driving functions
%   * Impulse responses
%   * 2.5D
%   * Tapering
% - Sound Field Simulations
% - Secondary Sources
% - Wave Field Synthesis (WFS)
%   * Pre-equalization
% - Spectral Division Method (SDM)
% - Near-Field Compensated Higher Order Ambisonics (NFC-HOA)
% - Local Sound Field Synthesis
% - Binaural Reproduction
%   * Headphone compensation
%   * SoundScape Renderer
% - Plotting
% - References
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


%% ===== Delayline =======================================================
% Delaying of time signals. This can be critical as very often the wanted delays
% are given as fractions of samples. This configuration section handles how
% those delays should be handled. As the default setting, integer only delays
% are used by rounding to the next larger integer delay.
% Beside choosing the actual delayline filter, the signal can also be resampled
% before delaying.
%
% Resample signal
%   'none'   - no resampling (default) 
%   'matlab' - use matlab's resample() function
%   'pm'     - use Parks-McClellan-Method to compute resample filter (firpm)
conf.delayline.resampling = 'none'; % / string
% Oversamplingfactor factor >= 1
% This should be in the order of (1/stepsize of fractional delays)
conf.delayline.resamplingfactor = 100; % / 1
% Order of Parks-McClellan resample filter (only for 'pm')
% This results in a filter length of resamplingfactor*resamplingorder
conf.delayline.resamplingorder = 64;
%
% Delayline filter
%   'integer'       - round to nearest integer delay (default)
%   'zoh'           - round to next larger integer delay
%   'lagrange'      - lagrange interpolator (FIR Filter)
%   'least_squares' - least squares FIR interpolation filter
%   'thiran'        - Thiran's allpass IIR filter
%   'farrow'        - use the Farrow structure (to be implemented)
conf.delayline.filter = 'integer';  % string
% Order of delayline filter (only for Lagrange, Least-Squares & Thiran)
conf.delayline.filterorder = 0;  % / 1
% Number of parallel filters in Farrow structure
% (only for 'farrow');
conf.delayline.filternumber = 1; % / 1


%% ===== Sound Field Synthesis (SFS) =====================================
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
%
% === Time Domain Implementation ===
% Adjust the starting time in time domain driving functions.
% This can be set to
%   'system'   - the first secondary source will be active at t=0
%   'source'   - the virtual source will be active at t=0
% Setting it to 'system' is most convenient when simulating single sources as
% you will always see activity in the sound field for t>0. Setting it to
% 'source' helps you to simulate different sources as you can time align them
% easily. Note, that for virtual sources outside of the array this can mean you
% will see no activity inside the listening area until the time has passed, that
% the virtual source needs from its position until the nearest secondary source.
% Also note, using 'source' for systems with unbounded listening areas -- e.g.
% linear arrays -- focussed virtual sources may not be placed arbitrarily
% far from the secondary sources.
conf.t0 = 'system'; % string
% Bandpass filter applied in sound_field_imp()
conf.usebandpass = false; % boolean
conf.bandpassflow = 10; % / Hz
conf.bandpassfhigh = 20000; % / Hz


%% ===== Sound Field Simulations =========================================
% Simulations of monochromatic or time domain sound field
%
% xyz-resolution for sound field simulations, this value is applied along every
% desired dimension, except if only one point is desired
conf.resolution = 300; % / samples
% Phase of omega of sound field (change this value to create monochromatic sound
% fields with different phases, for example this can be useful to create a movie)
conf.phase = 0; % / rad


% ===== Secondary Sources ================================================
% Settings of the used loudspeaker array
%
% Number of secondary sources
conf.secondary_sources.number = 64; % integer
% Diameter/Length of secondary source array
conf.secondary_sources.size = 3; % / m
% Center of array, X0
conf.secondary_sources.center = [0 0 0]; % / m
% Array geometry
% Possible values are: 'line', 'box', 'rounded-box', 'circle', 'sphere', 'custom'
conf.secondary_sources.geometry = 'circle'; % string
% exclusive for 'rounded-box' array geometry. Defines the bending radius for
% the corners of the smoothed box
conf.secondary_sources.corner_radius = 0.0; % / m
% Vector containing custom secondary source positions and directions.
% This is used if geometry = 'custom' is specified.
% conf.secondary_sources.x0 = [x0; y0; z0; nx0; ny0; nz0; weight];
% Or it could also be a SOFA struct or file name, in this case the positions are
% extracted from the provided SOFA file.
conf.secondary_sources.x0 = []; % / m
% Grid for the spherical array. Note, that you have to download and install the
% spherical grids from an additional source. For available grids see:
% http://github.com/sfstoolbox/data/tree/master/spherical_grids
% An exception are Gauss grids, which are available via 'gauss' and will be
% calculated on the fly allowing very high number of secondary sources.
conf.secondary_sources.grid = 'equally_spaced_points'; % string


%% ===== Wave Field Synthesis (WFS) ======================================
% Settings for WFS, see Spors et al. (2008) for an introduction
%
% === Pre-Equalization ===
% WFS can be implemented very efficiently using a delay-line with different
% amplitudes and convolving the whole signal once with the so called
% pre-equalization filter [References]. If we have aliasing in our system we
% only want to use the pre-equalization filter until the aliasing frequency,
% because of the energy the aliasing is adding to the spectrum above this
% frequency (which means the frequency response over the aliasing frequency is
% already "correct") [Reference]
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
% desired FIR filter order, results in N+1 taps
conf.wfs.hpreFIRorder = 128; % even integer


%% ===== Spectral Division Method (SDM) ==================================
% Settings for SDM, see Ahrens, Spors (2010) for an introduction
%
% Use the evanescent part of the driving function for SDM
conf.sdm.withev = true; % boolean


%% ===== Near-Field Compensated Higher Order Ambisonics (NFC-HOA) ========
% Settings for NFC-HOA, see Ahrens (2012) for an introduction
%
% Normally the order of NFC-HOA is set by the nfchoa_order() function which
% returns the highest order for which no aliasing occurs. If you wish to use
% another order you can set it manually here, otherwise leave it blank
conf.nfchoa.order = []; % integer
% Additional weighting of the modal coefficients by window function
conf.nfchoa.modal_window = 'rect';  % string
% Window type. Available windows are:
%   'rect'                     - all coefficients are weighted by 1.0
%   'kaiser', 'kaiser-bessel'  - Kaiser aka. Kaiser-Bessel window
conf.nfchoa.modal_window_parameter = 0.0;  % float
% Scalar parameter for window, if applicable. Effect for distinct window:
%   'rect'    - no effect
%   'kaiser'  - [0,inf]. trade-off between main-lobe width and side-lobe levels.
%               0.0 results in the rectangular window and the smallest main-lobe
%               width. infinity results in a dirac impulse.


%% ===== Local Sound Field Synthesis =====================================
% Settings for Local SFS, 
%
% Method the virtual secondary sources should be driven
conf.localsfs.method = 'wfs'; % 'wfs' or 'nfchoa'
conf.localsfs.usetapwin = false; % boolean
conf.localsfs.tapwinlen = 0.5; % 0..1
% WFS settings
conf.localsfs.wfs = conf.wfs;

% === Local Sound Field Synthesis using Virtual Secondary Sources ===
% see Spors, Ahrens (2010) for an introduction

% Virtual secondary sources (vss)
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0, 0];
conf.localsfs.vss.geometry = 'circular';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.grid = 'equally_spaced_points';
% driving function to create the focused sources, i.e. virtual secondary
% sources
conf.localsfs.vss.driving_functions = 'default';
%
% linear vss distribution: rotate the distribution orthogonal to the progation
% direction of the desired sound source
% circular vss distribution: truncate the distribution to a circular arc
% which satisfies the secondary source selection criterions ( source normal
% aligns with propagation directions of desired sound source )
conf.localsfs.vss.consider_target_field = true;
%
% vss distribution is further truncated if parts of it cannot be correctly
% reproduced, because they lie outside the area which is surrounded by the real
% loudspeakers (secondary sources)
conf.localsfs.vss.consider_secondary_sources = true;

% === Local Sound Field Synthesis using Spatial Bandwidth Limitation ===
conf.localsfs.sbl.order = 27;
conf.localsfs.sbl.fc = [];
conf.localsfs.sbl.Npw = [];
% driving function to create the plane waves
conf.localsfs.sbl.driving_functions = 'default';

%% ===== Binaural reproduction ===========================================
% Settings regarding all the stuff with impulse responses from the SFS_ir and
% SFS_binaural_synthesis folders
%
% Use interpolation to get the desired HRTF or BRIR for binaural simulation. If this
% is disabled, the HRTF/BRIR returned by a nearest neighbour search is used instead.
conf.ir.useinterpolation = true; % boolean
% You can choose the way the points for interpolation are selected. Depending on the
% geometry of the measured HRTF/BRIR data set, the interpolation will be done between
% two or three HRTFs. Available methods:
%   'nearestneighbour'  - Interpolation between nearest neighbours. This only works
%                         for interpolation points on a circle in the horizontal plane.
%   'delaunay'          - Interpolation between surrounding points according to
%                         Delaunay triangulation. This only works for interpolation
%                         points on a sphere. See validation script
%                         test_interpolation_point_selection.m for examples.
conf.ir.interpolationpointselection = 'nearestneighbour';
% You can choose between the following interpolation methods:
%   'simple'      - Interpolation in the time domain performed samplewise. This
%                   does not heed the times of arrival of the impulse responses.
%   'freqdomain'  - Interpolation in the frequency domain performed separately
%                   for magnitude and phase.
%                   This method cannot work properly if there is too much noise in
%                   the phase information at low frequencies which is often the
%                   case for measured HRTFs. Low frequencies can be corrected
%                   according to theory, see e.g. the corrected KEMAR HRTFs published
%                   at http://github.com/spatialaudio/lf-corrected-kemar-hrtfs.
%                   The implementation of this method suffers from circular shifting,
%                   see test_interpolation_methods.m in the validation folder. For
%                   typical HRIRs with leading and trailing zeros, the error is
%                   negligible.
conf.ir.interpolationmethod = 'simple';
%
% If you have HRIRs in the form of the SimpleFreeFieldHRIR SOFA convention, zeros
% are padded at the beginning of every impulse response corresponding to their
% measurement distance. If you know that your measured HRIRs already have a
% given pre-delay, add the pre-delay here and accordingly less zero padding will
% be applied. In this case you can lose samples from the beginning of the
% impulse response. If you are not sure, choose a value of 0.
conf.ir.hrirpredelay = 0; % / samples
%
% === Headphone compensation ===
% Headphone compensation
conf.ir.usehcomp = false; % boolean
% Headphone compensation file for left and right ear.
conf.ir.hcompfile = 'data/headphone_compensation/QU_KEMAR_AKGK601_hcomp.wav'; % string
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
% Normalize the sound field for plotting
conf.plot.usenormalisation = true; % boolean
% Normalisation method. Available methods are:
%   'auto'      - 'center' if center of sound field > 0.3, otherwise 'max'
%   'center'    - center of sound field == 1
%   'max'       - max of sound field == 1
conf.plot.normalisation = 'auto'; % string
% Plot mode (uses the GraphDefaults function). Available modes are:
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
% Matlab/Octave. For example to get the Matlab default colormap set 'jet'.
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


%% ===== References ======================================================
%
% Spors, Rabenstein, Ahrens - The Theory of Wave Field Synthesis Revisited, 124
% AES Convention, Paper 7358, 2008. http://bit.ly/ZCvyQ6
%
% Ahrens, Spors - Sound Field Reproduction Using Planar and Linear Arrays of
% Loudspeakers, Transactions on Audio, Speech, and Language Processing, p.
% 2038-50, 2010. http://bit.ly/10dpA9r
%
% Spors, Ahrens - Local Sound Field Synthesis by Virtual Secondary Sources, 40
% AES Conference, Paper 6-3, 2010. http://bit.ly/1t3842v
%
% Ahrens - Analytic Methods of Sound Field Synthesis. Springer, 2012.
