function boolean = test_localwfs()
%TEST_LOCALWFS tests the local WFS driving functions
%
%   Usage: boolean = test_localwfs()
%
%
%   Output parameters:
%       booelan - true or false
%
%   TEST_LOCALWFS() checks if the local WFS driving functions for the 
%   monochromatic case are working. Different sound fields are calculated and
%   plotted for visual inspection.

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


%% ===== Configuration ===================================================
boolean = false;
% Parameters
conf = SFS_config_example;
conf.plot.useplot = false;
conf.showprogress = true;
conf.resolution = 1000;
conf.plot.loudspeakers = true;
conf.plot.realloudspeakers = false;
conf.usetapwin = true;


%% ===== Circular secondary sources ======================================
% config for real loudspeaker array
conf.dimension = '2D';
conf.secondary_sources.geometry = 'circular';
conf.secondary_sources.number = 56;
conf.secondary_sources.size = 2;
conf.secondary_sources.center = [0, 0, 0];
conf.driving_functions = 'default';
conf.xref = conf.secondary_sources.center;
% listening area
X = [0 0 0];
xs = [1.0, -1.0, 0];  % propagation direction of plane wave
src = 'pw';
xrange = [-1, 1];
yrange = [-1, 1];
zrange = 0;
f = 7000;
% --- Circular virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0, 0];
conf.localsfs.vss.geometry = 'circular';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)
% --- Linear virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0.2, 0];
conf.localsfs.vss.geometry = 'linear';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)
% --- Box shaped virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0, 0];
conf.localsfs.vss.geometry = 'box';
conf.localsfs.vss.number = 4*56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)


%% =====  Linear secondary sources =======================================
% config for real loudspeaker array
conf.dimension = '2D';
conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.number = 56;
conf.secondary_sources.size = 2;
conf.secondary_sources.center = [0, 1, 0];
conf.driving_functions = 'default';
% listening area
X = [0 0 0];
conf.xref = X;
xs = [1.0, -1.0, 0];  % propagation direction of plane wave
src = 'pw';
xrange = [-1, 1];
yrange = [-1, 1];
zrange = 0;
f = 7000;
% --- Circular virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0, 0];
conf.localsfs.vss.geometry = 'circular';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)
% --- Linear virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0.2, 0];
conf.localsfs.vss.geometry = 'linear';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)
% --- Box shaped virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0, 0];
conf.localsfs.vss.geometry = 'box';
conf.localsfs.vss.number = 4*56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)


%% =====  Box shaped secondary sources ===================================
% config for real loudspeaker array
conf.dimension = '2D';
conf.secondary_sources.geometry = 'box';
conf.secondary_sources.number = 4*56;
conf.secondary_sources.size = 2;
conf.secondary_sources.center = [0, 0, 0];
conf.driving_functions = 'default';
% listening area
X = [0 0 0];
conf.xref = X;
xs = [1.0, -1.0, 0];  % propagation direction of plane wave
src = 'pw';
xrange = [-1, 1];
yrange = [-1, 1];
zrange = 0;
f = 7000;
% --- Circular virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0, 0];
conf.localsfs.vss.geometry = 'circular';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)
% --- Linear virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0.2, 0];
conf.localsfs.vss.geometry = 'linear';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)
% --- Box shaped virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0, 0];
conf.localsfs.vss.geometry = 'box';
conf.localsfs.vss.number = 4*56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)

boolean = true;
