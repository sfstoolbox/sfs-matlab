%% initialize
close all;
clear variables;
% SFS Toolbox
SFS_start;

%% Parameters
conf = SFS_config_example;

conf.dimension = '2.5D';
conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.number = 60;
conf.secondary_sources.size = 4;
conf.secondary_sources.center = [0, 0, 0];

conf.plot.useplot = false;

conf.scattering.Nse = 23;
conf.scattering.Nce = 23;

conf.showprogress = true;
conf.resolution = 400;

ns = [0, -1, 0];  % propagation direction of plane wave
xs = [0,  1.0, 0];  % position of point source
f = 1200;
xrange = [-2 2];
yrange = [-2 2];
zrange = 0;

% scatterer
sigma = inf;  % admittance of scatterer (inf to soft scatterer)
R = 0.3;
xq = [ 0, -1, 0];
xt = [ 0.5, 0.5, 0];
conf.xref = xq;
     
display(conf.scattering)

%% Spherical Expansion
% spherical expansion
A1sph = sphexp_mono_pw(ns,f,xq,conf);
A1sph_shift = sphexp_mono_pw(ns,f,xq+xt,conf);
A2sph = sphexp_mono_ps(xs,'R',f,xq,conf);
A2sph_shift = sphexp_mono_ps(xs,'R', f,xq+xt,conf);
% regular-to-regular spherical reexpansion (translatory shift)
[RRsph, RRsphm] = sphexp_mono_translation(-xt, 'RR', f, conf);
A1spht = RRsph*A1sph;
A2spht = RRsph*A2sph;

% Evaluate spherical basis functions 
[Jsphn, Hsphn, Ysphnm] = ...
  sphbasis_mono_XYZgrid(xrange,yrange,zrange,f,xq,conf);
%
[Jsphnt, Hsphnt, Ysphnmt] = ...
  sphbasis_mono_XYZgrid(xrange,yrange,zrange,f,xq+xt,conf);

% compute fields
P1sph = sound_field_mono_basis(A1sph, Jsphn, Ysphnm, conf);
P1sph_shift = sound_field_mono_basis(A1sph_shift, Jsphnt, Ysphnmt, conf);
P1spht = sound_field_mono_basis(A1spht, Jsphnt, Ysphnmt, conf);
P2sph = sound_field_mono_basis(A2sph, Jsphn, Ysphnm, conf);
P2sph_shift = sound_field_mono_basis(A2sph_shift, Jsphnt, Ysphnmt, conf);
P2spht = sound_field_mono_basis(A2spht, Jsphnt, Ysphnmt, conf);

% plot
[~,~,~,x1,y1,z1] = xyz_grid(xrange,yrange,zrange,conf);

plot_sound_field(P1sph ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('plane wave');
plot_sound_field(P1sph_shift ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('plane wave (shifted expansion)');
plot_sound_field(P1spht ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('plane wave (shifted reexpansion)');
plot_sound_field(P2sph ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('point source');
plot_sound_field(P2sph_shift ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('point source (shifted expansion)');
plot_sound_field(P2spht ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('point source (shifted reexpansion)');

%% Scattering with Sphere
% scatterer is assumed to be located at coordinates origin
B1sph = sphexp_mono_scatter(A1sph, R, sigma, f, conf); 
B2sph = sphexp_mono_scatter(A2sph, R, sigma, f, conf);

% compute fields
P1sph_scatter = sound_field_mono_basis(A1sph, Hsphn, Ysphnm, conf);
P2sph_scatter = sound_field_mono_basis(A2sph, Hsphn, Ysphnm, conf);

%% WFS Reproduction of Spherical Expansion
% regular spherical expansion of plane wave and point source
A1sph = sphexp_mono_pw(ns,f,xq,conf);
A2sph = sphexp_mono_ps(xs,'R',f,xq,conf);
% loudspeakers
x0 = secondary_source_positions(conf);
x0pw = secondary_source_selection(x0,ns,'pw');
x0ps = secondary_source_selection(x0,xs,'ps');
% compute driving functions and sound fields
D1sph = driving_function_mono_wfs_sphexp(x0pw(:,1:3),x0pw(:,4:6),A1sph,'R',f,xq,conf);
D2sph = driving_function_mono_wfs_sphexp(x0ps(:,1:3),x0ps(:,4:6),A2sph,'R',f,xq,conf);
% compute fields
P1sphwfs = sound_field_mono(xrange,yrange,zrange,x0pw,'ps',D1sph,f,conf);
P2sphwfs = sound_field_mono(xrange,yrange,zrange,x0ps,'ps',D2sph,f,conf);
% plot
[~,~,~,x1,y1,z1] = xyz_grid(xrange,yrange,zrange,conf);

plot_sound_field(P1sphwfs ,x1,y1,z1, x0pw, conf);
title('plane wave');
plot_sound_field(P2sphwfs ,x1,y1,z1, x0ps, conf);
title('point source');

%% WFS Reproduction of focused source using time reversal
% singular spherical expansion of point source
B1sph = sphexp_mono_ps(xs, 'S', f, xs, conf);
% loudspeakers
x0 = secondary_source_positions(conf);
x0ps = secondary_source_selection(x0,xs,'ps');
x0fs = secondary_source_selection(x0,[xs,0 1 0],'fs');
% compute driving functions
D1sph = driving_function_mono_wfs_sphexp(x0ps(:,1:3),x0ps(:,4:6),B1sph,'S',f,xs,conf);
D2sph = driving_function_mono_wfs_sphexp(x0fs(:,1:3),x0fs(:,4:6),B1sph,'S',f,xs,conf);
% compute fields
P1sphwfs = sound_field_mono(xrange,yrange,zrange,x0,'ps',D1sph,f,conf);
P2sphwfs = sound_field_mono(xrange,yrange,zrange,x0,'ps',conj(D2sph),f,conf);
% plot
[~,~,~,x1,y1,z1] = xyz_grid(xrange,yrange,zrange,conf);

plot_sound_field(P1sphwfs,x1,y1,z1, x0, conf);
title('point source');
plot_sound_field(P2sphwfs,x1,y1,z1, x0, conf);
title('focused source');

%% WFS Reproduction using virtual Scatterer and time reversal
% regular spherical expansion of plane wave and point source
A1sph = sphexp_mono_pw(ns,f,xq,conf);
A2sph = sphexp_mono_ps(xs,'R',f,xq,conf);
% compute timereversed incident field
A1sph_timereversed = sphexp_mono_timereverse(A1sph);
A2sph_timereversed = sphexp_mono_timereverse(A2sph);
% compute scattered field
B1sph = sphexp_mono_scatter(A1sph_timereversed, R, sigma, f, conf); 
B2sph = sphexp_mono_scatter(A2sph_timereversed, R, sigma, f, conf);
% loudspeakers
x0 = secondary_source_positions(conf);
x0pw = secondary_source_selection(x0,ns,'pw');
x0ps = secondary_source_selection(x0,xs,'ps');
% compute driving functions
D1sph = driving_function_mono_wfs_sphexp(x0pw(:,1:3),x0pw(:,4:6),B1sph,'S',f,xq,conf);
D2sph = driving_function_mono_wfs_sphexp(x0ps(:,1:3),x0ps(:,4:6),B2sph,'S',f,xq,conf);
% compute fields
P1sphwfs = sound_field_mono(xrange,yrange,zrange,x0pw,'ps',conj(D1sph),f,conf);
P2sphwfs = sound_field_mono(xrange,yrange,zrange,x0ps,'ps',conj(D2sph),f,conf);
% plot
[~,~,~,x1,y1,z1] = xyz_grid(xrange,yrange,zrange,conf);

plot_sound_field(P1sphwfs ,x1,y1,z1, x0pw, conf);
title('plane wave');
plot_scatterer(xq,R);
plot_sound_field(P2sphwfs ,x1,y1,z1, x0ps, conf);
title('point source');
plot_scatterer(xq,R);

%% generic NFCHOA
% loudspeakers
x0 = secondary_source_positions(conf);

A1spht_shift = RRsphm*sphexp_bandlimit(A1sph_shift,10);
A2spht_shift = RRsphm*sphexp_bandlimit(A2sph_shift,10);

% compute driving functions
D1sphhoa = driving_function_mono_nfchoa_sphexp(x0(:,1:3), A1sph, f, conf);
D2sphhoa = driving_function_mono_nfchoa_sphexp(x0(:,1:3), A2sph, f, conf);
% compute driving functions for shifted expansions
D1sphhoat = driving_function_mono_nfchoa_sphexp(x0(:,1:3), A1spht_shift, f, conf);
D2sphhoat = driving_function_mono_nfchoa_sphexp(x0(:,1:3), A2spht_shift, f, conf);

% compute fields
P1sphhoa = sound_field_mono(xrange,yrange,zrange,x0,'ps',D1sphhoa,f,conf);
P2sphhoa = sound_field_mono(xrange,yrange,zrange,x0,'ps',D2sphhoa,f,conf);
% compute fields
P1sphhoat = sound_field_mono(xrange,yrange,zrange,x0,'ps',D1sphhoat,f,conf);
P2sphhoat = sound_field_mono(xrange,yrange,zrange,x0,'ps',D2sphhoat,f,conf);

% plot
[~,~,~,x1,y1,z1] = xyz_grid(xrange,yrange,zrange,conf);

plot_sound_field(P1sphhoa, x1,y1,z1, x0, conf);
title('NFCHOA: plane wave');
plot_sound_field(P2sphhoa ,x1,y1,z1, x0, conf);
title('NFCHOA: point source');
plot_sound_field(P1sphhoat/5, x1,y1,z1, x0, conf);
title('NFCHOA: plane wave (shifted reexpansion)');
plot_sound_field(P2sphhoat/5, x1,y1,z1, x0, conf);
title('NFCHOA: point source (shifted reexpansion)');