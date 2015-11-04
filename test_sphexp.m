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
conf.resolution = 100;

ns = [0, -1, 0];  % propagation direction of plane wave
xs = [0,  2.0, 0];  % position of point source

X = [-2 2];
Y = [-2 2];
Z = 0;

f = 500;
Nse = 10;

% scatterer
sigma = inf;  % admittance of scatterer (inf to soft scatterer)
R = 0.3;
xq = [ 0, -1, 0];
xt = [ 0.5, 0.5, 0];
conf.xref = xq;

%% Spherical Expansion
% spherical expansion
A1sph = sphexp_mono_pw(ns,Nse,f,xq,conf);
A1sph_shift = sphexp_mono_pw(ns,Nse,f,xq+xt,conf);
A2sph = sphexp_mono_ps(xs,'R',Nse, f,xq,conf);
A2sph_shift = sphexp_mono_ps(xs,'R', Nse, f,xq+xt,conf);
% regular-to-regular spherical reexpansion (translatory shift)
[RRsph, RRsphm] = sphexp_mono_translation(xt, 'RR', Nse, f, conf);
A1spht = RRsph*A1sph;
A2spht = RRsph*A2sph;

% Evaluate spherical basis functions 
[jn, h2n, Ynm] = ...
  sphbasis_mono_grid(X,Y,Z,Nse,f,xq,conf);
%
[jnt, ~, Ynmt] = ...
  sphbasis_mono_grid(X, Y, Z,Nse,f,xq+xt,conf);

% compute fields
P1sph = sound_field_mono_sphbasis(A1sph, jn, Ynm);
P1sph_shift = sound_field_mono_sphbasis(A1sph_shift, jnt, Ynmt);
P1spht = sound_field_mono_sphbasis(A1spht, jnt, Ynmt);
P2sph = sound_field_mono_sphbasis(A2sph, jn, Ynm);
P2sph_shift = sound_field_mono_sphbasis(A2sph_shift, jnt, Ynmt);
P2spht = sound_field_mono_sphbasis(A2spht, jnt, Ynmt);

plot_sound_field(P1sph ,X, Y, Z, [], conf);
plot_scatterer(xq,R);
title('plane wave');
plot_sound_field(P1sph_shift ,X, Y, Z, [], conf);
plot_scatterer(xq,R);
title('plane wave (shifted expansion)');
plot_sound_field(P1spht ,X, Y, Z, [], conf);
plot_scatterer(xq,R);
title('plane wave (shifted reexpansion)');
plot_sound_field(P2sph ,X, Y, Z, [], conf);
plot_scatterer(xq,R);
title('point source');
plot_sound_field(P2sph_shift ,X, Y, Z, [], conf);
plot_scatterer(xq,R);
title('point source (shifted expansion)');
plot_sound_field(P2spht ,X, Y, Z, [], conf);
plot_scatterer(xq,R);
title('point source (shifted reexpansion)');

%% Scattering with Sphere
% scatterer is assumed to be located at coordinates origin
B1sph = sphexp_mono_scatter(A1sph, R, sigma, f, conf); 
B2sph = sphexp_mono_scatter(A2sph, R, sigma, f, conf);

% compute fields
P1sph_scatter = sound_field_mono_sphbasis(B1sph, h2n, Ynm);
P2sph_scatter = sound_field_mono_sphbasis(B2sph, h2n, Ynm);