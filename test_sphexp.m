%% initialize
close all;
clear variables;
% SFS Toolbox
SFS_start;

%% Parameters
conf = SFS_config;

conf.dimension = '2.5D';

conf.plot.useplot = false;

conf.showprogress = true;
conf.resolution = 100;

ns = [0, -1, 0];  % propagation direction of plane wave
xs = [0,  2.0, 0];  % position of point source

X = [-2 2];
Y = [-2 2];
Z = 0;

f = 1000;
Nse = 20;

% scatterer
sigma = inf;  % admittance of scatterer (inf to soft scatterer)
R = 0.3;
xq = [ 0, 0, 0];
xt = [ 0.5, 0, 0];
conf.xref = xq;


%% Spherical Expansion
% regular spherical expansion at xq
Apwsph = sphexp_mono_pw(ns,Nse,f,xq,conf);
Apssph = sphexp_mono_ps(xs,'R',Nse, f,xq,conf);
% regular spherical expansion at xq+xt (x' = x - xt)
Apwsph_t = sphexp_mono_pw(ns,Nse,f,xq+xt,conf);
Apssph_t = sphexp_mono_ps(xs,'R', Nse, f,xq+xt,conf);
% regular-to-regular spherical reexpansion (translatory shift of -xt)
[RRsph, RRsphm] = sphexp_mono_translation(-xt, 'RR', Nse, f, conf);
A1sph_re = RRsph*Apwsph_t;
A2sph_re = RRsph*Apssph_t;

% Evaluate spherical basis functions 
[jn, h2n, Ynm] = ...
  sphbasis_mono_grid(X,Y,Z,Nse,f,xq,conf);
%
[jnt, ~, Ynmt] = ...
  sphbasis_mono_grid(X, Y, Z,Nse,f,xq+xt,conf);

% compute fields
Ppwsph = sound_field_mono_sphbasis(Apwsph, jn, Ynm);
Ppwsph_t = sound_field_mono_sphbasis(Apwsph_t, jnt, Ynmt);
Ppwsph_re = sound_field_mono_sphbasis(A1sph_re, jn, Ynm);
Ppssph = sound_field_mono_sphbasis(Apssph, jn, Ynm);
Ppssph_t = sound_field_mono_sphbasis(Apssph_t, jnt, Ynmt);
Ppssph_re = sound_field_mono_sphbasis(A2sph_re, jn, Ynm);

plot_sound_field(Ppwsph ,X, Y, Z, [], conf);
plot_scatterer(xq,R);
title('plane wave');
plot_sound_field(Ppwsph_t ,X, Y, Z, [], conf);
plot_scatterer(xq,R);
title('plane wave (shifted expansion)');
plot_sound_field(Ppwsph_re ,X, Y, Z, [], conf);
plot_scatterer(xq,R);
title('plane wave (shifted reexpansion)');
plot_sound_field(Ppssph ,X, Y, Z, [], conf);
plot_scatterer(xq,R);
title('point source');
plot_sound_field(Ppssph_t ,X, Y, Z, [], conf);
plot_scatterer(xq,R);
title('point source (shifted expansion)');
plot_sound_field(Ppssph_re ,X, Y, Z, [], conf);
plot_scatterer(xq,R);
title('point source (shifted reexpansion)');

%% Scattering with Sphere
% scatterer is assumed to be located at coordinates origin
B1sph = sphexp_mono_scatter(Apwsph, R, sigma, f, conf); 
B2sph = sphexp_mono_scatter(Apssph, R, sigma, f, conf);

% compute fields
Ppwsph_scatter = sound_field_mono_sphbasis(B1sph, h2n, Ynm);
Ppssph_scatter = sound_field_mono_sphbasis(B2sph, h2n, Ynm);