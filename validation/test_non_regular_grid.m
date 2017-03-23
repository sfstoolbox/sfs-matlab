function status = test_non_regular_grid(modus)
%TEST_IMPULSE_RESPONSES tests time behavior of WFS and local WFS
%
%   Usage: status = test_impulse_responses(modus)
%
%   Input parameters:
%       modus   - 0: numerical
%                 1: visual
%
%   Output parameters:
%       status  - true or false
%
%   TEST_IMPULSE_RESPONSES(modus) compares the time-frequency response of
%   WFS and local WFS by calculating impulse responses, their frequency
%   spectrum, and spatial-temporal sound field.

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
%% Parameters
conf = SFS_config;
conf.showprogress = true;
conf.resolution = 400;
if modus
    conf.plot.useplot = true;
    conf.plot.loudspeakers = true;
    conf.plot.realloudspeakers = false;
    conf.plot.usedb = false;
end
conf.tapwinlen = 0.3;
% config for array
conf.dimension = '2.5D';
conf.secondary_sources.geometry = 'circular';
conf.secondary_sources.number = 56;
conf.secondary_sources.size = 3;
conf.secondary_sources.center = [0, 0, 0];
conf.driving_functions = 'default';
conf.xref = [0,0,0];
% listening area, virtual source
xs = [0.0, 2.5, 0];  % propagation direction of plane wave
src = 'ps';
f = 1000;
t = 190/conf.fs;

conf.usenormalisation = true;

%% ===== Computation =====================================================
% regular grid
Xreg = [-1.5 1.5];
Yreg = [-1, 1.55];
Zreg = 0;

% non regular grid
alpha = 2*pi / 360 * (0:360-1);
r = linspace(0,conf.secondary_sources.size/2,50);
[alpha,r] = ndgrid(alpha,r);

Xnon  = r.*cos(alpha);
Ynon  = r.*sin(alpha);
Znon = 0;

% sound fields
conf.plot.normalisation = 'center';
[~] = sound_field_mono_wfs(Xreg,Yreg,Zreg,xs,src,f,conf);
[~] = sound_field_mono_wfs(Xnon,Ynon,Znon,xs,src,f,conf);

conf.plot.normalisation = 'max';
[~] = sound_field_imp_wfs(Xreg,Yreg,Zreg,xs,src,t,conf);
[~] = sound_field_imp_wfs(Xnon,Ynon,Znon,xs,src,t,conf);


status = true;
