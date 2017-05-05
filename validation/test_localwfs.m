function status = test_localwfs(modus)
%TEST_LOCALWFS tests the local WFS driving functions
%
%   Usage: status = test_localwfs(modus)
%
%   Input parameters:
%       modus   - 0: numerical
%                 1: visual
%
%   Output parameters:
%       status  - true or false
%
%   TEST_LOCALWFS(modus) checks if the local WFS driving functions for the
%   monochromatic case are working. Different sound fields are calculated and
%   plotted for visual inspection.

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


status = false;


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);


%% ===== Configuration ===================================================
% Parameters
conf = SFS_config;
if modus
    conf.plot.useplot = true;
    conf.plot.normalisation = 'center';
    conf.plot.loudspeakers = true;
    conf.plot.realloudspeakers = false;
end
conf.showprogress = true;
conf.resolution = 1000;
conf.usetapwin = true;


%% ===== Circular secondary sources ======================================
% config for real loudspeaker array
conf.dimension = '2D';
conf.secondary_sources.geometry = 'circular';
conf.secondary_sources.number = 56;
conf.secondary_sources.size = 2;
conf.secondary_sources.center = [0 0 0];
conf.driving_functions = 'default';
conf.xref = conf.secondary_sources.center;
% listening area
X = [0 0 0];
xs = [1.0 -1.0 0];  % propagation direction of plane wave
src = 'pw';
X = [-1 1];
Y = [-1 1];
Z = 0;
f = 7000;
% --- Circular virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0 0 0];
conf.localsfs.vss.geometry = 'circular';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
[~] = sound_field_mono_localwfs_vss(X,Y,Z,xs,src,f,conf);
% --- Linear virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0 0.2 0];
conf.localsfs.vss.geometry = 'linear';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
[~] = sound_field_mono_localwfs_vss(X,Y,Z,xs,src,f,conf);
% --- Box shaped virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0 0 0];
conf.localsfs.vss.geometry = 'box';
conf.localsfs.vss.number = 4*56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
[~] = sound_field_mono_localwfs_vss(X,Y,Z,xs,src,f,conf);


%% =====  Linear secondary sources =======================================
% config for real loudspeaker array
conf.dimension = '2D';
conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.number = 56;
conf.secondary_sources.size = 2;
conf.secondary_sources.center = [0 1 0];
conf.driving_functions = 'default';
% listening area
X = [0 0 0];
conf.xref = X;
xs = [1.0 -1.0 0];  % propagation direction of plane wave
src = 'pw';
X = [-1 1];
Y = [-1 1];
Z = 0;
f = 7000;
% --- Circular virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0 0 0];
conf.localsfs.vss.geometry = 'circular';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
[~] = sound_field_mono_localwfs_vss(X, Y, Z,xs,src,f,conf);
% --- Linear virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0 0.2 0];
conf.localsfs.vss.geometry = 'linear';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
[~] = sound_field_mono_localwfs_vss(X, Y, Z,xs,src,f,conf);
% --- Box shaped virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0 0 0];
conf.localsfs.vss.geometry = 'box';
conf.localsfs.vss.number = 4*56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
[~] = sound_field_mono_localwfs_vss(X, Y, Z,xs,src,f,conf);


%% =====  Box shaped secondary sources ===================================
% config for real loudspeaker array
conf.dimension = '2D';
conf.secondary_sources.geometry = 'box';
conf.secondary_sources.number = 4*56;
conf.secondary_sources.size = 2;
conf.secondary_sources.center = [0 0 0];
conf.driving_functions = 'default';
% listening area
X = [0 0 0];
conf.xref = X;
xs = [1.0 -1.0 0];  % propagation direction of plane wave
src = 'pw';
X = [-1 1];
Y = [-1 1];
Z = 0;
f = 7000;
% --- Circular virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0 0 0];
conf.localsfs.vss.geometry = 'circular';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
[~] = sound_field_mono_localwfs_vss(X,Y,Z,xs,src,f,conf);
% --- Linear virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0 0.2 0];
conf.localsfs.vss.geometry = 'linear';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
[~] = sound_field_mono_localwfs_vss(X,Y,Z,xs,src,f,conf);
% --- Box shaped virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0 0 0];
conf.localsfs.vss.geometry = 'box';
conf.localsfs.vss.number = 4*56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
[~] = sound_field_mono_localwfs_vss(X,Y,Z,xs,src,f,conf);


status = true;
