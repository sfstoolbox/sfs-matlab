function boolean = test_non_regular_grid()
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
conf.plot.usedb = false;
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
tau = 190;

conf.usenormalisation = true;

%% ===== Computation =====================================================
% regular grid
Xreg = [-1.5 1.5];
Yreg = [-1, 1.55];
Zreg = 0;

% non regular grid
alpha = 2*pi / 360 * (0:360-1);
r = linspace(0, conf.secondary_sources.size/2, 50);
[alpha, r] = ndgrid(alpha,r);

Xnon  = r.*cos(alpha);
Ynon  = r.*sin(alpha);
Znon = 0;

% sound fields
conf.plot.normalisation = 'center';
sound_field_mono_wfs(Xreg,Yreg,Zreg, xs, src, f, conf);
sound_field_mono_wfs(Xnon,Ynon,Znon, xs, src, f, conf);

conf.plot.normalisation = 'max';
sound_field_imp_wfs(Xreg,Yreg,Zreg, xs, src, tau, conf);
sound_field_imp_wfs(Xnon,Ynon,Znon, xs, src, tau, conf);

boolean = true;
