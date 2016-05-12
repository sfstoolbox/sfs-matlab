function boolean = test_impulse_responses()
%TEST_IMPULSE_RESPONSES tests time behavior of WFS and local WFS
%
%   Usage: boolean = test_impulse_responses()
%
%   Output parameters:
%       booelan - true or false
%
%   TEST_IMPULSE_RESPONSES() compares the time-frequency response of
%   WFS and local WFS by calculating impulse responses, their frequency
%   spectrum, and spatial-temporal sound field.

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Team                                   *
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


%% ===== Configuration ===================================================
boolean = false;
%% Parameters
conf = SFS_config;
conf.showprogress = true;
conf.resolution = 400;
conf.plot.useplot = true;
conf.plot.loudspeakers = true;
conf.plot.realloudspeakers = false;
conf.plot.usedb = true;
conf.tapwinlen = 0.3;
% config for virtual array
conf.localsfs.method = 'wfs';
conf.localsfs.wfs = conf.wfs;
conf.localsfs.usetapwin = true;
conf.localsfs.tapwinlen = 0.3;
conf.localsfs.vss.size = 0.6;
conf.localsfs.vss.center = [-1, 0, 0];
conf.localsfs.vss.geometry = 'circular';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
% config for real array
conf.dimension = '2.5D';
conf.secondary_sources.geometry = 'circular';
conf.secondary_sources.number = 56;
conf.secondary_sources.size = 3;
conf.secondary_sources.center = [0, 0, 0];
conf.driving_functions = 'default';
conf.xref = conf.localsfs.vss.center;
% impulse response
conf.ir.usehcomp = false;
% listening area, virtual source
xs = [0.0, 2.5, 0];  % propagation direction of plane wave
src = 'ps';
X = [-1.55 1.55];
Y = [-1.55, 1.55];
Z = 0;


%% ===== Computation =====================================================
%% temporal impulse responses
irs = dummy_irs(conf.N,conf);

% === WFS ===
% calculate impulse response
s_wfs = ir_wfs(conf.xref,pi/2,xs,src,irs,conf);
% plot frequency response
[S_wfs, ~, f_wfs] = easyfft(s_wfs(:,1)./max(abs(s_wfs(:,1))), conf);
% plot spatio-temporal sound field
sound_field_imp_wfs(X,Y,Z, xs, src, 190, conf);

% === Local WFS ===
conf.tapwinlen = 1.0;
% calculate prefilter
[conf.wfs.hpreflow, conf.wfs.hprefhigh] = ...
  localwfs_findhpref(conf.xref, pi/2, xs, src, conf);
conf.localsfs.wfs = conf.wfs;
% calculate impulse response
s_lwfs = ir_localwfs(conf.xref,pi/2,xs,src,irs,conf);
% plot frequency response
[S_lwfs, ~, f_lwfs] = easyfft(s_lwfs(:,1)./max(abs(s_lwfs(:,1))), conf);
% plot spatio-temporal sound field
sound_field_imp_localwfs(X,Y,Z, xs, src, 360, conf);

boolean = true;
