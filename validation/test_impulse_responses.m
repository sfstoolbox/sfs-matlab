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
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
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
%% Parameters
conf = SFS_config_example;
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

% just to test the new delayline implementation
conf.fracdelay.pre.method = 'none';
conf.fracdelay.pre.resample.method = 'matlab';
conf.fracdelay.pre.resample.factor = 50;
conf.fracdelay.filter = 'lagrange';
conf.fracdelay.length = 1;

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

% % === Local WFS ===
% conf.tapwinlen = 1.0;
% % calculate prefilter
% [conf.wfs.hpreflow, conf.wfs.hprefhigh] = ...
%   localwfs_findhpref(conf.xref, pi/2, xs, src, conf);
% conf.localsfs.wfs = conf.wfs;
% % calculate impulse response
% s_lwfs = ir_localwfs(conf.xref,pi/2,xs,src,irs,conf);
% % plot frequency response
% [S_lwfs, ~, f_lwfs] = easyfft(s_lwfs(:,1)./max(abs(s_lwfs(:,1))), conf);
% % plot spatio-temporal sound field
% sound_field_imp_localwfs(X,Y,Z, xs, src, 360, conf);

boolean = true;
