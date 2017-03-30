function status = test_plotting(modus)
%TEST_PLOT tests the correctness of test_plot()
%
%   Usage: status = test_plot(modus)
%
%   Input parameters:
%       modus   - 0: numerical (not available)
%                 1: visual
%
%   Output parameters:
%       status  - true or false
%
%   TEST_PLOT(modus) creates plots for monochromatic and time-domain sound field
%   with different dimensions and different grids.

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
t = 300/conf.fs;


%% ===== Monochromatic plots =============================================
conf.plot.normalisation = 'center';
% 3-D plots
tmp = conf.resolution;
conf.resolution = 50;
sound_field_mono_wfs([-2 2],[-2 2],[-2 2],xs,src,f,conf)
title('3D monochromatic sound field')
conf.resolution = tmp;
% 2-D plots
sound_field_mono_wfs([-2 2],[-2 2],0,xs,src,f,conf)
title('2D monochromatic sound field, xy-axes')
sound_field_mono_wfs([-2 2],0,[-2 2],xs,src,f,conf)
title('2D monochromatic sound field, xz-axes')
sound_field_mono_wfs(0,[-2 2],[-2 2],xs,src,f,conf)
title('2D monochromatic sound field, yz-axes')
% 1-D plots
sound_field_mono_wfs([-2 2],0,0,xs,src,f,conf)
title('1D monochromatic sound field, x-axis')
sound_field_mono_wfs(0,[-2 2],0,xs,src,f,conf)
title('1D monochromatic sound field, y-axis')
sound_field_mono_wfs(0,0,[-2 2],xs,src,f,conf)
title('1D monochromatic sound field, z-axis')
% === Custom grid ===
x1 = randi([-2000 2000],125000,1)/1000;
x2 = randi([-2000 2000],125000,1)/1000;
x3 = randi([-2000 2000],125000,1)/1000;
% 3-D plots
sound_field_mono_wfs(x1,x2,x3,xs,src,f,conf)
title('3D monochromatic sound field, custom grid')
% 2-D plots
sound_field_mono_wfs(x1,x2,0,xs,src,f,conf)
title('2D monochromatic sound field, xy-axes, custom grid')
sound_field_mono_wfs(x1,0,x2,xs,src,f,conf)
title('2D monochromatic sound field, xz-axes, custom grid')
sound_field_mono_wfs(0,x1,x2,xs,src,f,conf)
title('2D monochromatic sound field, yz-axes, custom grid')
% 1-D plots
sound_field_mono_wfs(x1,0,0,xs,src,f,conf)
title('1D monochromatic sound field, x-axis, custom grid')
sound_field_mono_wfs(0,x2,0,xs,src,f,conf)
title('1D monochromatic sound field, y-axis, custom grid')
sound_field_mono_wfs(0,0,x3,xs,src,f,conf)
title('1D monochromatic sound field, z-axis, custom grid')


%% ===== Time-domain plots ===============================================
conf.plot.normalisation = 'max';
% 3-D plots
tmp = conf.resolution;
conf.resolution = 50;
sound_field_imp_wfs([-2 2],[-2 2],[-2 2],xs,src,t,conf)
title('3D time-domain sound field')
conf.resolution = tmp;
% 2-D plots
sound_field_imp_wfs([-2 2],[-2 2],0,xs,src,t,conf)
title('2D time-domain sound field, xy-axes')
sound_field_imp_wfs([-2 2],0,[-2 2],xs,src,t,conf)
title('2D time-domain sound field, xz-axes')
sound_field_imp_wfs(0,[-2 2],[-2 2],xs,src,t,conf)
title('2D time-domain sound field, yz-axes')
% 1-D plots
sound_field_imp_wfs([-2 2],0,0,xs,src,t,conf)
title('1D time-domain sound field, x-axis')
sound_field_imp_wfs(0,[-2 2],0,xs,src,t,conf)
title('1D time-domain sound field, y-axis')
sound_field_imp_wfs(0,0,[-2 2],xs,src,t,conf)
title('1D time-domain sound field, z-axis')
% === Custom grid ===
x1 = randi([-2000 2000],125000,1)/1000;
x2 = randi([-2000 2000],125000,1)/1000;
x3 = randi([-2000 2000],125000,1)/1000;
% 3-D plots
sound_field_imp_wfs(x1,x2,x3,xs,src,t,conf)
title('3D time-domain sound field, custom grid')
% 2-D plots
sound_field_imp_wfs(x1,x2,0,xs,src,t,conf)
title('2D time-domain sound field, xy-axes, custom grid')
sound_field_imp_wfs(x1,0,x2,xs,src,t,conf)
title('2D time-domain sound field, xz-axes, custom grid')
sound_field_imp_wfs(0,x1,x2,xs,src,t,conf)
title('2D time-domain sound field, yz-axes, custom grid')
% 1-D plots
sound_field_imp_wfs(x1,0,0,xs,src,t,conf)
title('1D time-domain sound field, x-axis, custom grid')
sound_field_imp_wfs(0,x2,0,xs,src,t,conf)
title('1D time-domain sound field, y-axis, custom grid')
sound_field_imp_wfs(0,0,x3,xs,src,t,conf)
title('1D time-domain sound field, z-axis, custom grid')


%% ===== Time-domain plots in dB =========================================
conf.plot.usedb = true;
% 3-D plots
tmp = conf.resolution;
conf.resolution = 50;
sound_field_imp_wfs([-2 2],[-2 2],[-2 2],xs,src,t,conf)
title('3D time-domain sound field')
conf.resolution = tmp;
% 2-D plots
sound_field_imp_wfs([-2 2],[-2 2],0,xs,src,t,conf)
title('2D time-domain sound field, xy-axes')
sound_field_imp_wfs([-2 2],0,[-2 2],xs,src,t,conf)
title('2D time-domain sound field, xz-axes')
sound_field_imp_wfs(0,[-2 2],[-2 2],xs,src,t,conf)
title('2D time-domain sound field, yz-axes')
% 1-D plots
sound_field_imp_wfs([-2 2],0,0,xs,src,t,conf)
title('1D time-domain sound field, x-axis')
sound_field_imp_wfs(0,[-2 2],0,xs,src,t,conf)
title('1D time-domain sound field, y-axis')
sound_field_imp_wfs(0,0,[-2 2],xs,src,t,conf)
title('1D time-domain sound field, z-axis')
% === Custom grid ===
x1 = randi([-2000 2000],125000,1)/1000;
x2 = randi([-2000 2000],125000,1)/1000;
x3 = randi([-2000 2000],125000,1)/1000;
% 3-D plots
sound_field_imp_wfs(x1,x2,x3,xs,src,t,conf)
title('3D time-domain sound field, custom grid')
% 2-D plots
sound_field_imp_wfs(x1,x2,0,xs,src,t,conf)
title('2D time-domain sound field, xy-axes, custom grid')
sound_field_imp_wfs(x1,0,x2,xs,src,t,conf)
title('2D time-domain sound field, xz-axes, custom grid')
sound_field_imp_wfs(0,x1,x2,xs,src,t,conf)
title('2D time-domain sound field, yz-axes, custom grid')
% 1-D plots
sound_field_imp_wfs(x1,0,0,xs,src,t,conf)
title('1D time-domain sound field, x-axis, custom grid')
sound_field_imp_wfs(0,x2,0,xs,src,t,conf)
title('1D time-domain sound field, y-axis, custom grid')
sound_field_imp_wfs(0,0,x3,xs,src,t,conf)
title('1D time-domain sound field, z-axis, custom grid')


status = true;
