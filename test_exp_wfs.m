%% initialize
close all;
clear variables;
% SFS Toolbox
SFS_start;

%% Parameters
conf = SFS_config_example;
conf.showprogress = true;
conf.dimension = '2D';

% plotting
conf.plot.usedb = false;
conf.plot.useplot = false;
conf.usenormalisation = true;
conf.resolution = 300;

X = [-2 2];
Y = [-4 0];
Z = 0;

% secondary sources
conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.number = 64;
conf.secondary_sources.size = 4;
conf.secondary_sources.center = [0, 0, 0];

f = 1000;
ns = [0, -1, 0];  % propagation direction of plane wave
xs = [0, 2, 0];  % position of point source
ls = xs;          % position of line source

xq = [0 -2 0];  % expansion center
conf.xref = xq;  % reference position

xt = [ 0.5, 0.5, 0];  % position of the sweet spot
rt = 1.0;  % "size" of the sweet spot
Nce = circexp_truncation_order(rt, f, 1e-6, conf);  % 2.5DHOA order for sweet spot

sigma = inf;  % (inf for sound soft scatterer)

%% Spherical Expansion Coefficients
% regular spherical expansion at xq
Apwnm = sphexp_mono_pw(ns, Nce,f,xq,conf);
Apsnm = sphexp_mono_ps(xs,'R', Nce,f,xq,conf);
Alsnm = sphexp_mono_ls(xs,'R', Nce,f,xq,conf);
% regular spherical expansion at xq+xt
Apwnm_t = sphexp_mono_pw(ns,Nce,f,xq+xt,conf);
Apsnm_t = sphexp_mono_ps(xs,'R',Nce,f,xq+xt,conf);
Alsnm_t = sphexp_mono_ls(xs,'R',Nce,f,xq+xt,conf);
% compute timereversed incident field
Apwnm_rev = sphexp_mono_timereverse(Apwnm);
Apsnm_rev = sphexp_mono_timereverse(Apsnm);
Alsnm_rev = sphexp_mono_timereverse(Alsnm);
% compute timereversed incident field
Apwnm_t_rev = sphexp_mono_timereverse(Apwnm_t);
Apsnm_t_rev = sphexp_mono_timereverse(Apsnm_t);
Alsnm_t_rev = sphexp_mono_timereverse(Alsnm_t);
% scatterer is assumed to be located at the expansion center (i.e. xq)
Bpwnm = sphexp_mono_scatter(Apwnm_rev, rt, sigma, f, conf); 
Bpsnm = sphexp_mono_scatter(Apsnm_rev, rt, sigma, f, conf);
Blsnm = sphexp_mono_scatter(Alsnm_rev, rt, sigma, f, conf);
% scatterer is assumed to be located at the expansion center (i.e. xq+xt)
Bpwnm_t = sphexp_mono_scatter(Apwnm_t_rev, rt, sigma, f, conf); 
Bpsnm_t = sphexp_mono_scatter(Apsnm_t_rev, rt, sigma, f, conf);
Blsnm_t = sphexp_mono_scatter(Alsnm_t_rev, rt, sigma, f, conf);

%% Circular Expansion Coefficients
% regular circular expansion at xq
Apwm = circexp_mono_pw(ns, Nce,f,xq,conf);
Alsm = circexp_mono_ls(xs, 'R', Nce,f,xq,conf);
% regular circular expansion at xq+xt
Apwm_t = circexp_mono_pw(ns,Nce,f,xq+xt,conf);
Alsm_t = circexp_mono_ls(xs, 'R', Nce,f,xq+xt,conf);
% compute timereversed incident field
Apwm_rev = circexp_mono_timereverse(Apwm);
Alsm_rev = circexp_mono_timereverse(Alsm);
% compute timereversed incident field
Apwm_t_rev = circexp_mono_timereverse(Apwm_t);
Alsm_t_rev = circexp_mono_timereverse(Alsm_t);
% scatterer is assumed to be located at the expansion center (i.e. xq)
Bpwm = circexp_mono_scatter(Apwm_rev, rt, sigma, f, conf); 
Blsm = circexp_mono_scatter(Alsm_rev, rt, sigma, f, conf);
% scatterer is assumed to be located at the expansion center (i.e. xq+xt)
Bpwm_t = circexp_mono_scatter(Apwm_t_rev, rt, sigma, f, conf); 
Blsm_t = circexp_mono_scatter(Alsm_t_rev, rt, sigma, f, conf);

%% generic WFS in spatial domain using Spherical Expansion Coefficients
conf.dimension = '2D';
% loudspeakers (TODO: implicit selection of loudspeakers in driving
% function)
x0 = secondary_source_positions(conf);
% compute driving functions
Dpw = driving_function_mono_wfs_sphexp(x0(:,1:3),x0(:,4:6), Apwnm,'R',f,xq,conf);
Dps = driving_function_mono_wfs_sphexp(x0(:,1:3),x0(:,4:6), Apsnm,'R',f,xq,conf);
Dls = driving_function_mono_wfs_sphexp(x0(:,1:3),x0(:,4:6), Alsnm,'R',f,xq,conf);
% compute fields
Ppw = sound_field_mono(X,Y,Z,x0,'ps',Dpw,f,conf);
Pps = sound_field_mono(X,Y,Z,x0,'ps',Dps,f,conf);
Pls = sound_field_mono(X,Y,Z,x0,'ps',Dls,f,conf);
% plot
plot_sound_field(Ppw, X, Y, Z, [], conf);
title('2.5D WFS with spherical expansion (spatial domain): plane wave');
plot_sound_field(Pps, X, Y, Z, [], conf);
title('2.5D WFS with spherical expansion (spatial domain): point source');
plot_sound_field(Pls, X, Y, Z, [], conf);
title('2.5D WFS with spherical expansion (spatial domain): line source');

%% generic WFS in spatial domain using circular Expansion Coefficients
conf.dimension = '2D';
% loudspeakers (TODO: implicit selection of loudspeakers in driving
% function)
x0 = secondary_source_positions(conf);
% compute driving functions
Dpw = driving_function_mono_wfs_circexp(x0(:,1:3),x0(:,4:6), Apwm,'R',f,xq,conf);
Dls = driving_function_mono_wfs_circexp(x0(:,1:3),x0(:,4:6), Alsm,'R',f,xq,conf);
% compute fields
Ppw = sound_field_mono(X,Y,Z,x0,'ls',Dpw,f,conf);
Pls = sound_field_mono(X,Y,Z,x0,'ls',Dls,f,conf);
% plot
plot_sound_field(Ppw, X, Y, Z, [], conf);
title('2D WFS with circular expansion (spatial domain): plane wave');
plot_sound_field(Pls, X, Y, Z, [], conf);
title('2D WFS with circular expansion (spatial domain): line source');

%% LWFS in spatial domain using virtual Spherical Scatterer and time reversal
conf.dimension = '2D';
% loudspeakers (TODO: implicit selection of loudspeakers in driving
% function)
x0 = secondary_source_positions(conf);
% compute driving functions
Dpw = driving_function_mono_wfs_sphexp(x0(:,1:3),x0(:,4:6), Bpwnm,'S',f,xq,conf);
Dps = driving_function_mono_wfs_sphexp(x0(:,1:3),x0(:,4:6), Bpsnm,'S',f,xq,conf);
Dls = driving_function_mono_wfs_sphexp(x0(:,1:3),x0(:,4:6), Blsnm,'S',f,xq,conf);
% compute fields
Ppw = sound_field_mono(X,Y,Z,x0,'ps',conj(Dpw),f,conf);
Pps = sound_field_mono(X,Y,Z,x0,'ps',conj(Dps),f,conf);
Pls = sound_field_mono(X,Y,Z,x0,'ps',conj(Dls),f,conf);
% plot
plot_sound_field(Ppw, X, Y, Z, [], conf);
title('2.5D local WFS with spherical expansion (spatial domain): plane wave');
plot_sound_field(Pps, X, Y, Z, [], conf);
title('2.5D local WFS with spherical expansion (spatial domain): point source');
plot_sound_field(Pls, X, Y, Z, [], conf);
title('2.5D local WFS with spherical expansion (spatial domain): line source');

%% LWFS in spatial domain using virtual Spherical Scatterer and time reversal
conf.dimension = '2D';
% loudspeakers (TODO: implicit selection of loudspeakers in driving
% function)
x0 = secondary_source_positions(conf);
% compute driving functions
Dpw = driving_function_mono_wfs_circexp(x0(:,1:3),x0(:,4:6),Bpwm,'S',f,xq,conf);
Dls = driving_function_mono_wfs_circexp(x0(:,1:3),x0(:,4:6),Blsm,'S',f,xq,conf);
% compute fields
Ppw = sound_field_mono(X,Y,Z,x0,'ls',conj(Dpw),f,conf);
Pls = sound_field_mono(X,Y,Z,x0,'ls',conj(Dls),f,conf);
% plot
plot_sound_field(Ppw, X, Y, Z, [], conf);
title('2D local WFS with circular expansion (spatial domain): plane wave');
plot_sound_field(Pls, X, Y, Z, [], conf);
title('2D local WFS with circular expansion (spatial domain): line source');