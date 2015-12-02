function boolean = test_colormaps()
%TEST_COLORMAPS does sound field plots with different colormaps
%
%   Usage: boolean = test_colormaps()
%
%   Output parameters:
%       booelan - true or false
%
%   TEST_COLORMAPS() creates plots for monochromatic and time-domain sound field in dB
%   with different colormaps.

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


%% ===== Checking of input  parameters ===================================
nargmin = 0;
nargmax = 0;
narginchk(nargmin,nargmax);
boolean = false;

%% ===== Configuration ===================================================
conf = SFS_config_example;
xs = [0 -1 0];
src = 'pw';
f = 1000;
t = 300;
X = [-2 2];
Y = [-2 2];
Z = 0;

%% ===== Monochromatic plots =============================================
conf.plot.normalisation = 'center';
conf.plot.usedb = true;

color_maps = { ...
    'chromajs'; ...
    'gray'; ...
    };
color_maps_reversed = { ...
    'magma'; ...
    'viridis'; ...
    'cubehelix'; ...
    'inferno'; ...
    'bone'; ...
    };

[P,~,~,~,x0] = sound_field_mono_wfs([-2 2],[-2 2],0,xs,src,f,conf);
p = sound_field_imp_wfs(X,Y,Z,xs,src,t,conf);

for ii=1:length(color_maps)
    conf.plot.colormap = color_maps{ii};
    plot_sound_field(P,X,Y,Z,x0,conf)
    title(sprintf('Monochromatic, %s',color_maps{ii}))
    plot_sound_field(p,X,Y,Z,x0,conf)
    title(sprintf('Time-domain, %s',color_maps{ii}))
end
for ii=1:length(color_maps_reversed)
    conf.plot.colormap = color_maps_reversed{ii};
    plot_sound_field(P,X,Y,Z,x0,conf)
    colormap(flipud(colormap))
    title(sprintf('Monochromatic, %s reversed',color_maps_reversed{ii}))
    plot_sound_field(p,X,Y,Z,x0,conf)
    colormap(flipud(colormap))
    title(sprintf('Time-domain, %s reversed',color_maps_reversed{ii}))
end

boolean = true;
