%% initialize
close all;
clear variables;
% SFS Toolbox
SFS_start;

%% Parameters
conf = SFS_config_example;

conf.dimension = '2.5D';

conf.plot.useplot = false;

conf.showprogress = true;
conf.resolution = 200;

ns = [0, -1, 0];  % propagation direction of plane wave
xs = [0,  3.0, 0];  % position of line source

xrange = [-2 2];
yrange = [-2 2];
zrange = 0;

f = 1000;

xq = [ 0, 0, 0];

conf.xref = xq;

xt = [ 0.5, 0, 0];  % position of the sweet spot
rt = 1.0;  % "size" of the sweet spot
Nce = circexp_truncation_order(rt, f, 1e-6, conf);  % 2.5DHOA order for sweet spot
Ncet = circexp_truncation_order(norm(xt)+rt, f, 1e-6, conf);  % 2.5DHOA order for translation

%% Circular Expansion Coefficients
% regular circular expansion at xq
Apwm = circexp_mono_pw(ns, Nce,f,xq,conf);
Alsm = circexp_mono_ls(xs, 'R', Nce,f,xq,conf);
% regular circular expansion at xq+xt
Apwm_t = circexp_mono_pw(ns,Nce,f,xq+xt,conf);
Alsm_t = circexp_mono_ls(xs, 'R', Nce,f,xq+xt,conf);
% regular-to-regular cylindrical reexpansion (translatory shift)
[RR, RRm] = circexp_mono_translation(xt,'RR',Nce,f,conf);
Apwm_re = RRm*Apwm;
Alsm_re = RRm*Alsm;

%% Sound Fields
% Evaluate spherical basis functions on regular grid
[Jm, H2m, Ym] = ...
  circbasis_mono_grid(xrange,yrange,zrange,Nce,f,xq,conf);
%
[Jm_t, H2m_t, Ym_t] = ...
  circbasis_mono_grid(xrange,yrange,zrange,Nce,f,xq+xt,conf);
% compute fields for expansion at xq
Ppwm = sound_field_mono_circbasis(Apwm, Jm, Ym);
Plsm = sound_field_mono_circbasis(Alsm, Jm, Ym);
% compute fields for expansion at xq + xt
Ppwm_t = sound_field_mono_circbasis(Apwm_t, Jm_t, Ym_t);
Plsm_t = sound_field_mono_circbasis(Alsm_t, Jm_t, Ym_t);
% compute field for cylindrical reexpansion (translatory shift)
Ppwm_re = sound_field_mono_circbasis(Apwm_re, Jm, Ym_t);
Plsm_re = sound_field_mono_circbasis(Alsm_re, Jm, Ym_t);
% plot
[~,~,~,x1,y1,z1] = xyz_grid(xrange,yrange,zrange,conf);

plot_sound_field(Ppwm ,x1,y1,z1, [], conf);
plot_scatterer(xq,rt);
title('plane wave');
plot_sound_field(Ppwm_t ,x1,y1,z1, [], conf);
plot_scatterer(xq,rt);
title('plane wave (shifted expansion)');
plot_sound_field(Ppwm_re ,x1,y1,z1, [], conf);
plot_scatterer(xq,rt);
title('plane wave (re-expansion)');
plot_sound_field(Plsm ,x1,y1,z1, [], conf);
plot_scatterer(xq,rt);
title('line source');
plot_sound_field(Plsm_t,x1,y1,z1, [], conf);
plot_scatterer(xq,rt);
title('line source (shifted expansion)');
plot_sound_field(Plsm_re ,x1,y1,z1, [], conf);
plot_scatterer(xq,rt);
title('line source (reexpansion)');