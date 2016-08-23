function status = test_colormaps(modus)
%TEST_COLORMAPS does sound field plots with different colormaps
%
%   Usage: status = test_colormaps(modus)
%
%   Input parameters:
%       modus   - 0: numerical (not available)
%                 1: visual
%
%   Output parameters:
%       status  - true or false
%
%   TEST_COLORMAPS(modus) creates plots for monochromatic and time-domain sound
%   field in dB with different colormaps.

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


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
if ~modus
    warning('%s: numerical modus not available.',upper(mfilename));
    status = true;
    return;
end


%% ===== Configuration ===================================================
conf = SFS_config;
xs = [0 -1 0];
src = 'pw';
f = 1000;
t = 0.007;
X = [-2 2];
Y = [-2 2];
Z = 0;


%% ===== Sound field plots ===============================================
conf.plot.normalisation = 'center';
conf.plot.usedb = true;

color_maps = { ...
    'yellowred'; ...
    'gray'; ...
    };
color_maps_reversed = { ...
    'magma'; ...
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


status = true;
