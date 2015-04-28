%% initialize
close all;
clear variables;
% SFS Toolbox
SFS_start;

%% Parameters
conf = SFS_config_example;

conf.dimension = '2.5D';
conf.secondary_sources.geometry = 'circular';
conf.secondary_sources.number = 60;
conf.secondary_sources.size = 3;
conf.secondary_sources.center = [0, 0, 0];

conf.plot.useplot = false;

conf.scattering.Nse = 23;
conf.scattering.Nce = 23;

conf.showprogress = true;
conf.resolution = 400;

ns = [0, -1, 0];  % propagation direction of plane wave
xs = [0,  2, 0];  % position of point source
f = 1200;
xrange = [-2 2];
yrange = [-2 2];
zrange = 0;

xq = conf.secondary_sources.center;
xt = [ 0.5, 0.5, 0];
conf.xref = xq;
     
display(conf.scattering)

%% generic NFCHOA
% regular spherical expansion of plane wave and point source
A1sph = sphexp_mono_pw(ns,f,xq+xt,conf);
A2sph = sphexp_mono_ps(xs,'R', f,xq+xt,conf);
% regular-to-regular spherical reexpansion (translatory shift)
[RRsph, RRsphm] = sphexp_mono_translation(-xt, 'RR', f, conf);
% shift spherical expansion back to xq
A1sph_shift = RRsph*sphexp_bandlimit(A1sph,10);
A2sph_shift = RRsph*sphexp_bandlimit(A2sph,10);
% loudspeakers
x0 = secondary_source_positions(conf);
% compute driving functions
D1sphhoa = driving_function_mono_nfchoa_sphexp(x0(:,1:3), A1sph, f, conf);
D2sphhoa = driving_function_mono_nfchoa_sphexp(x0(:,1:3), A2sph, f, conf);
% compute driving functions for shifted expansions
D1sphhoa_shift = driving_function_mono_nfchoa_sphexp(x0(:,1:3), A1sph_shift, f, conf);
D2sphhoa_shift = driving_function_mono_nfchoa_sphexp(x0(:,1:3), A2sph_shift, f, conf);
% compute fields
P1sphhoa = sound_field_mono(xrange,yrange,zrange,x0,'ps',D1sphhoa,f,conf);
P2sphhoa = sound_field_mono(xrange,yrange,zrange,x0,'ps',D2sphhoa,f,conf);
% compute fields
P1sphhoa_shift = sound_field_mono(xrange,yrange,zrange,x0,'ps',D1sphhoa_shift,f,conf);
P2sphhoa_shift = sound_field_mono(xrange,yrange,zrange,x0,'ps',D2sphhoa_shift,f,conf);
% plot
[~,~,~,x1,y1,z1] = xyz_grid(xrange,yrange,zrange,conf);

plot_sound_field(P1sphhoa, x1,y1,z1, x0, conf);
title('NFCHOA: plane wave');
plot_sound_field(P2sphhoa ,x1,y1,z1, x0, conf);
title('NFCHOA: point source');
plot_sound_field(P1sphhoa_shift/5, x1,y1,z1, x0, conf);
title('NFCHOA: plane wave (shifted reexpansion)');
plot_sound_field(P2sphhoa_shift/5, x1,y1,z1, x0, conf);
title('NFCHOA: point source (shifted reexpansion)');