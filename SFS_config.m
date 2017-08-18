function conf = SFS_config()
%SFS_CONFIG returns a struct with the default configuration settings
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
%   Please don't edit this function to change your configuration settings!
%
%   See also: SFS_start

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


%% ===== Checking of input  parameters ===================================
nargmin = 0;
nargmax = 0;
narginchk(nargmin,nargmax);


%% ===== Table of Contents ===============================================
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
% - Local Sound Field Synthesis (LSFS)
% - Binaural Reproduction
%   * Headphone compensation
%   * SoundScape Renderer
% - Plotting
% - References


%% ===== Misc ============================================================
conf.tmpdir = '/tmp/sfs'; % string
% Debugging level. We are supporting two levels:
%   0 - normal mode
%   1 - debug modus, showing interim results and plots
conf.debug = 0; % 0 or 1
% Show a progress bar in selected loops (for example sound_field_mono). This can
% be useful, if you are using a high number of secondary sources.
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
% those delays should be handled. Beside choosing the actual delayline filter,
% the signal can also be resampled before delaying. As the default setting,
% integer only delays are used by rounding to the next larger integer delay.
% If you want to use a fractional delayline, a setting with a high accuracy is:
%   conf.delayline.resampling = 'pm';
%   conf.delayline.resamplingfactor = 8;
%   conf.delayline.resamplingorder = 64;
%   conf.delayline.filter = 'lagrange';
%   conf.delayline.filterorder = 9;
% Note, that the necessary interpolation accuracy highly depends on the
% actual use case and parametrisation, compare Winter, Spors (2016).
%
% Delayline filter
%   'integer'       - round to nearest integer delay (default)
%   'zoh'           - round to next larger integer delay
%   'lagrange'      - lagrange interpolator (FIR Filter)
%   'least_squares' - least squares FIR interpolation filter
%   'thiran'        - Thiran's allpass IIR filter
%   'farrow'        - use the Farrow structure (to be implemented)
conf.delayline.filter = 'integer';  % string
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
% functions are described at http://sfstoolbox.org. For possible other flags
% have a look into the driving functions as they can be quite specific.
% Most users can safely use the 'default' flag here.
conf.driving_functions = 'default'; % string
%
% === Impulse responses ===
% Length of impulse responses used in the time domain driving functions
% and for the creation of the binaural simulations.
% Don't worry, SFS checks for you if conf.N is large enough.
conf.N = 2048; % samples
%
% === 2.5D ===
% The amplitude will be correct at the point xref for 2.5D
% synthesis.
conf.xref = [0 0 0]; % / m
%
% === Tapering ===
% The truncation of the loudspeaker array leads to diffraction of the
% synthesized sound field. It has been shown that the truncation can be
% described by cylindrical waves originating from the edges of the array, see
% Sect. 8.3.2 in Born, Wolf (1999) for the general principle and Sect. 3.2 in
% Wierstorf (2014) for how it relates to WFS. Therefore a good method to reduce
% artifacts due to the diffraction edge waves is to fade out the amplitude of
% the driving function at the edges of the array. This method is called tapering
% and implemented using a Hann window in the SFS Toolbox.
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
%
% === Modal Weighting ===
% Additional weighting of the modal coefficients by window function used 
% for Near-Field Compensated Higher Order Ambisonics (NFC-HOA) and Local Sound
% Field Synthesis using Spatial Bandwidth Limitation (LSFS-SBL).
%
% Window type. Available windows are:
%   'rect'                     - all coefficients are weighted by 1.0
%   'max-rE'                   - 2D max-rE weighting
%   'kaiser', 'kaiser-bessel'  - Kaiser aka. Kaiser-Bessel window
%   'tukey'                    - modified Tukey (tapered cosine) window
conf.modal_window = 'rect';  % string
% Scalar parameter for window, if applicable. Effect for distinct window:
%   'rect'    - no effect
%   'max-rE'  - no effect
%   'kaiser'  - [0,inf]. trade-off between main-lobe width and side-lobe levels.
%               0.0 results in the rectangular window and the smallest main-lobe
%               width. infinity results in a dirac impulse.
%   'tukey'   - [0,1]. width of cosine tapering relative to the modal order. 
%               0.0 results in the rectangular window, 1.0 results in a modified
%               Hann-window
conf.modal_window_parameter = 0.0;  % float



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
% Exclusive for 'rounded-box' array geometry. Defines the bending radius for
% the corners of the smoothed box
conf.secondary_sources.corner_radius = 0.0; % / m
% Vector containing custom secondary source positions and directions.
% This is used if geometry = 'custom' is specified.
% conf.secondary_sources.x0 = [x0; y0; z0; nx0; ny0; nz0; weight];
% Or it could also be a SOFA struct or file name, in this case the positions are
% extracted from the provided SOFA file.
conf.secondary_sources.x0 = []; % / m
% Grid for a spherical array. Available grids are:
%   'equally_spaced_points' - Sphere with equal distance between grid points
%   'gauss'                 - Gauss grid
%   'fabian'                - grid of 3D HRTF measurement, available at
%                             http://dx.doi.org/10.14279/depositonce-5718
%
% Note, that 'equally_spaced_points' and 'fabian' are precomputed grids that
% will be automatically downloaded and cached on your disk. All available number
% of secondary sources for those grids can be seen at:
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
% pre-equalization filter, see Spors, Ahrens (2010a).
% Use WFS preequalization-filter
conf.wfs.usehpre = true; % boolean
% FIR or IIR pre-equalization filter
% NOTE: only FIR is working under octave at the moment
conf.wfs.hpretype = 'FIR'; % 'FIR' or 'IIR'
% Lower frequency limit of preequalization filter
% If we have a finite length (<10m) of the secondary source distribution we
% will have a 3dB increase at the very low frequencies and the pre-equalization
% filter should only start above those frequencies, see Sect. 7.2 in Spors,
% Ahrens (2010b).
conf.wfs.hpreflow = 50; % / Hz
% Upper frequency limit of preequalization filter
% If we have aliasing in our system we only want to use
% the pre-equalization filter until the aliasing frequency, because of the
% energy the aliasing is adding to the spectrum above this frequency, see
% Sect. 6.3 in Spors, Ahrens (2010a).
conf.wfs.hprefhigh = 1200; % / Hz
% IIR bandwidth for the Lagrange interpolation region
conf.wfs.hpreBandwidth_in_Oct = 2; % / octaves
% desired IIR filter order
conf.wfs.hpreIIRorder = 4; % integer
% desired FIR filter order, results in N+1 taps
conf.wfs.hpreFIRorder = 128; % even integer


%% ===== Spectral Division Method (SDM) ==================================
% Settings for SDM, see Ahrens, Spors (2010b) for an introduction
%
% Use the evanescent part of the driving function for SDM
conf.sdm.withev = true; % boolean


%% ===== Near-Field Compensated Higher Order Ambisonics (NFC-HOA) ========
% Settings for NFC-HOA, see Ahrens (2012) for an introduction
%
% Highest order used with NFC-HOA. If this is set to [], band-limited NFC-HOA is
% used and the order is set by nfchoa_order() which returns the highest order
% for which no aliasing occurs.
conf.nfchoa.order = []; % integer


%% ===== Local Wave Field Synthesis (LWFS) ===============================
% Settings for Local WFS
%
% === Local Wave Field Synthesis using Virtual Secondary Sources (LWFS-VSS)
% See Spors, Ahrens (2010) for an introduction.
%
% Method the virtual secondary sources should be driven
conf.localwfs_vss.method = 'wfs'; % 'wfs' or 'nfchoa'
% WFS settings for virtual secondary sources
conf.localwfs_vss.wfs = conf.wfs;
% Tapering of virtual secondary sources (only applied for WFS)
conf.localwfs_vss.usetapwin = false; % boolean
conf.localwfs_vss.tapwinlen = 0.5; % 0..1
% NFC-HOA settings for virtual secondary sources
conf.localwfs_vss.nfchoa = conf.nfchoa;
% Virtual secondary sources (see also: conf.secondary_sources)
conf.localwfs_vss.size = 0.4; % / m
conf.localwfs_vss.center = [0, 0, 0]; % / m
conf.localwfs_vss.geometry = 'circular'; % string
conf.localwfs_vss.number = 56; % integer
conf.localwfs_vss.grid = 'equally_spaced_points'; % string
% Driving function for virtual secondary sources
conf.localwfs_vss.driving_functions = 'default'; % string
% Linear VSS distribution: rotate the distribution orthogonal to the progation
% direction of the desired sound source.
% Circular VSS distribution: truncate the distribution to a circular arc
% which satisfies the secondary source selection criterions (source normal
% aligns with propagation directions of desired sound source).
conf.localwfs_vss.consider_target_field = true; % boolean
% VSS distribution is further truncated if parts of it cannot be correctly
% reproduced, because they lie outside the area which is surrounded by the real
% loudspeakers (secondary sources)
conf.localwfs_vss.consider_secondary_sources = true; % boolean
%
% === Local Wave Field Synthesis using Spatial Bandwidth Limitation (LWFS-SBL)
% See Hahn, Winter, Spors (2016) for an introduction.
% The centre of the local synthesis region is set by conf.xref
%
% Maximum modal order aka. spatial bandwidth of desired sound field. If left
% empty, the value is set by nfchoa_order(), which may suboptimal depending on
% the geometry, e.g. number of secondary sources and shape of the secondary 
% source distribution.
conf.localwfs_sbl.order = []; % integer
% Due to stability issues for the time-domain implementation of synthesis 
% of a point source, conventional WFS has to be used for the low frequencies. 
% fc defines the crossover frequency between the WFS and LWFS-SBL. If left
% empty, this frequency is estimated by aliasing_frequency().
conf.localwfs_sbl.fc = []; % float
% The spatially bandwidth-limited sound field is converted into plane wave
% decomposition which is then synthesised using conventional WFS for each
% individual plane wave. Npw defines the number of plane waves with their 
% directions distributed equi-angularly on the unit circle. If left empty,
% it is estimated based on the sampling frequency and size of the secondary 
% source distribution.
conf.localwfs_sbl.Npw = []; % integer


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
% Matlab/Octave. For example to get the Matlab default colormap set 'parula'.
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
% Ahrens (2012) - "Analytic Methods of Sound Field Synthesis", Springer,
% https://doi.org/10.1007/978-3-642-25743-8
%
% Ahrens, Spors (2010) - "Sound Field Reproduction Using Planar and Linear
% Arrays of Loudspeakers", Transactions on Audio, Speech, and Language
% Processing, vol. 18, no. 8, pp. 2038-2050,
% https://doi.org/10.1109/TASL.2010.2041106
%
% Born, Wolf (1999) - "Principles of Optics", Cambridge University Press, 7th
% edition.
%
% Hahn, Winter, Spors (2016) - "Local Wave Field Synthesis by Spatial
% Band-limitation in the Circular/Spherical Harmonics Domain", in 140th
% Convention of the Audio Engineering Society, Paper 9596,
% http://www.aes.org/e-lib/browse.cfm?elib=18294
%
% Spors, Ahrens (2010a) - "Analysis and Improvement of Pre-Equalization in
% 2.5-Dimensional Wave Field Synthesis", in 128th Convention of the Audio
% Engineering Society, Paper 8121,
% http://www.aes.org/e-lib/browse.cfm?elib=15418
%
% Spors, Ahrens (2010b) - "Local Sound Field Synthesis by Virtual Secondary
% Sources", in 40th Conference of the Audio Engineering Society, Paper 6-3,
% http://www.aes.org/e-lib/browse.cfm?elib=15561
%
% Spors, Rabenstein, Ahrens (2008) - "The Theory of Wave Field Synthesis
% Revisited", in 124th Convention of the Audio Engineering Society, Paper 7358,
% http://www.aes.org/e-lib/browse.cfm?elib=14488
%
% Wierstorf (2014) - "Perceptual Assessment of Sound Field Synthesis",
% TU Berlin, https://doi.org/10.14279/depositonce-4310
%
% Winter, Spors (2016) - "On fractional delay interpolation for local wave
% field synthesis", 24th European Signal Processing Conference (EUSIPCO),
% pp. 2415-2419, https://doi.org/10.1109/EUSIPCO.2016.7760682
