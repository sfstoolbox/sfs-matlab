%% initialize
% close all;
clear variables;
% SFS Toolbox
SFS_start;

%% Parameters
conf = SFS_config_example;
conf.showprogress = true;

% plotting
conf.plot.useplot = false;
conf.resolution = 400;

xrange = [-2 2];
yrange = [-2 2];
zrange = 0;

[~,~,~,x1,y1,z1] = xyz_grid(xrange,yrange,zrange,conf);

% secondary sources
conf.dimension = '2.5D';
conf.secondary_sources.geometry = 'circular';
conf.secondary_sources.number = 200;
conf.secondary_sources.size = 3;
conf.secondary_sources.center = [0, 0, 0];

ns = [0, -1, 0];  % propagation direction of plane wave
xs = [0,  2, 0];  % position of point source
f = 300;
Nse = 10;

xq = conf.secondary_sources.center;
xt = [ 0.5, 0.5, 0];
conf.xref = xq;
r0 = conf.secondary_sources.size / 2;

%% Spherical Expansion Coefficients
% regular spherical expansion of plane wave and point source at xq
A1nm_original = sphexp_mono_pw(ns,Nse,f,xq,conf);
A2nm_original = sphexp_mono_ps(xs,'R',Nse,f,xq,conf);
% regular spherical expansion of plane wave and point source at xq+xt
A1nm = sphexp_mono_pw(ns,Nse,f,xq+xt,conf);
A2nm = sphexp_mono_ps(xs,'R',Nse,f,xq+xt,conf);
% regular-to-regular spherical reexpansion (translatory shift)
% [RRsph, RRsphm] = sphexp_mono_translation(-xt, 'RR', f, conf);
% shift spherical expansion back to xq
% A1sph_shift = RRsph*sphexp_bandlimit(A1sph,10);
% A2sph_shift = RRsph*sphexp_bandlimit(A2sph,10);

%% generic NFCHOA in spatial domain
% loudspeakers
x0 = secondary_source_positions(conf);
% compute driving functions
D1 = driving_function_mono_nfchoa_sphexp(x0(:,1:3), A1nm_original, f, conf);
D2 = driving_function_mono_nfchoa_sphexp(x0(:,1:3), A2nm_original, f, conf);
% compute fields
P1 = sound_field_mono(xrange,yrange,zrange,x0,'ps',D1,f,conf);
P2 = sound_field_mono(xrange,yrange,zrange,x0,'ps',D2,f,conf);
% plot
plot_sound_field(P1, x1,y1,z1, [], conf);
title('NFCHOA (spatial domain): plane wave');
plot_sound_field(P2, x1,y1,z1, [], conf);
title('NFCHOA (spatial domain): point source');

%% generice NFCHOA in spherical harmonics domain
% compute driving functions
D1nm = driving_function_mono_nfchoa_sht_sphexp(A1nm_original,f,conf);
D2nm = driving_function_mono_nfchoa_sht_sphexp(A2nm_original,f,conf);
% compute fields
P1 = sound_field_mono_nfchoa_sht(xrange,yrange,zrange, D1nm, f, conf);
P2 = sound_field_mono_nfchoa_sht(xrange,yrange,zrange, D2nm, f, conf);
% plot
plot_sound_field(P1, x1,y1,z1, [], conf);
title('NFCHOA (sht domain): plane wave');
plot_sound_field(P2, x1,y1,z1, [], conf);
title('NFCHOA (sht domain): point source');