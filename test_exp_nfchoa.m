%% initialize
% close all;
clear variables;
% SFS Toolbox
SFS_start;

%% Parameters
conf = SFS_config_example;
conf.showprogress = true;

% plotting
conf.plot.usedb = false;
conf.plot.useplot = false;
conf.usenormalisation = false;
conf.resolution = 400;

xrange = [-2 2];
yrange = [-2 2];
zrange = 0;

[~,~,~,x1,y1,z1] = xyz_grid(xrange,yrange,zrange,conf);

% secondary sources
conf.secondary_sources.geometry = 'circular';
conf.secondary_sources.number = 56;
conf.secondary_sources.size = 3;
conf.secondary_sources.center = [0, 0, 0];

f = 1000;
ns = [0, -1, 0];  % propagation direction of plane wave
xs = [0,  2, 0];  % position of point source

xq = conf.secondary_sources.center;  % expansion center
conf.xref = xq;  % reference position

xt = [ 0.5, 0.5, 0];  % position of the sweet spot
rt = 0.3;  % "size" of the sweet spot
Ncet = circexp_truncation_order(norm(xt)+rt, f, 1e-6, conf);  % 2.5DHOA order for translation
Nce = circexp_truncation_order(rt, f, 1e-6, conf);  % 2.5DHOA order for sweet spot
Nse = sphexp_truncation_order(rt, f, 1e-6, conf);  % 3DHOA order for sweet spot 

%% Spherical Expansion Coefficients
% regular spherical expansion at xq
Apwnm_25D = sphexp_mono_pw(ns, Nce,f,xq,conf);
Apsnm_25D = sphexp_mono_ps(xs,'R', Nce,f,xq,conf);
Alsnm_25D = sphexp_mono_ls(xs,'R', Nce,f,xq,conf);
% regular spherical expansion at xq
Apwnm_3D = sphexp_mono_pw(ns, Nse,f,xq,conf);
Apsnm_3D = sphexp_mono_ps(xs,'R', Nse,f,xq,conf);
Alsnm_3D = sphexp_mono_ls(xs,'R', Nse,f,xq,conf);
% regular spherical expansion at xq+xt
Apwnm = sphexp_mono_pw(ns,Ncet,f,xq+xt,conf);
Apsnm = sphexp_mono_ps(xs,'R',Ncet,f,xq+xt,conf);
Alsnm = sphexp_mono_ls(xs,'R',Ncet,f,xq+xt,conf);
% regular-to-regular spherical reexpansion (translatory shift)
[RRsph, RRsphm] = sphexp_mono_translation(-xt, 'RR', Ncet, f, conf);
% shift spherical expansion back to xq
Apwnm_shift = RRsph*sphexp_truncate(Apwnm, Nce);
Apsnm_shift = RRsph*sphexp_truncate(Apsnm, Nce);
Alsnm_shift = RRsph*sphexp_truncate(Alsnm, Nce); 

%% generic NFCHOA in spatial domain
conf.dimension = '2.5D';
% loudspeakers
x0 = secondary_source_positions(conf);
% compute driving functions
Dpw = driving_function_mono_nfchoa_sphexp(x0(:,1:3), Apwnm_25D, f, conf);
Dps = driving_function_mono_nfchoa_sphexp(x0(:,1:3), Apsnm_25D, f, conf);
Dls = driving_function_mono_nfchoa_sphexp(x0(:,1:3), Alsnm_25D, f, conf);
% compute driving functions from shifted fields
Dpw_shift = driving_function_mono_nfchoa_sphexp(x0(:,1:3), Apwnm_shift, f, conf);
Dps_shift = driving_function_mono_nfchoa_sphexp(x0(:,1:3), Apsnm_shift, f, conf);
Dls_shift = driving_function_mono_nfchoa_sphexp(x0(:,1:3), Alsnm_shift, f, conf);
% compute fields
Ppw = sound_field_mono(xrange,yrange,zrange,x0,'ps',Dpw,f,conf);
Pps = sound_field_mono(xrange,yrange,zrange,x0,'ps',Dps,f,conf);
Pls = sound_field_mono(xrange,yrange,zrange,x0,'ls',Dls,f,conf);
% compute fields from shifted driving functions
Ppw_shift = sound_field_mono(xrange,yrange,zrange,x0,'ps',Dpw_shift,f,conf);
Pps_shift = sound_field_mono(xrange,yrange,zrange,x0,'ps',Dps_shift,f,conf);
Pls_shift = sound_field_mono(xrange,yrange,zrange,x0,'ls',Dls_shift,f,conf);
% plot
plot_sound_field(Ppw, x1,y1,z1, [], conf);
title('2.5D NFCHOA (spatial domain): plane wave');
plot_sound_field(Pps, x1,y1,z1, [], conf);
title('2.5D NFCHOA (spatial domain): point source');
plot_sound_field(Pls, x1,y1,z1, [], conf);
title('2.5D NFCHOA (spatial domain): line source');
plot_sound_field(Ppw_shift, x1,y1,z1, [], conf);
title('2.5D NFCHOA (spatial domain): plane wave (shifted expansion)');
plot_sound_field(Pps_shift, x1,y1,z1, [], conf);
title('2.5D NFCHOA (spatial domain): point source (shifted expansion)');
plot_sound_field(Pls_shift, x1,y1,z1, [], conf);
title('2.5D NFCHOA (spatial domain): line source (shifted expansion)');

%% generic 2.5D NFCHOA in spherical harmonics domain
conf.dimension = '2.5D';
% compute sht of driving functions
Dpwnm = driving_function_mono_nfchoa_sht_sphexp(Apwnm_25D,f,conf);
Dpsnm = driving_function_mono_nfchoa_sht_sphexp(Apsnm_25D,f,conf);
Dlsnm = driving_function_mono_nfchoa_sht_sphexp(Alsnm_25D,f,conf);
% compute sht of driving functions from shifted fields
Dpwnm_shift = driving_function_mono_nfchoa_sht_sphexp(Apwnm_shift, f, conf);
Dpsnm_shift = driving_function_mono_nfchoa_sht_sphexp(Apsnm_shift, f, conf);
Dlsnm_shift = driving_function_mono_nfchoa_sht_sphexp(Alsnm_shift, f, conf);
% compute fields
Ppw = sound_field_mono_nfchoa_sht(xrange,yrange,zrange, Dpwnm, f, conf);
Pps = sound_field_mono_nfchoa_sht(xrange,yrange,zrange, Dpsnm, f, conf);
Pls = sound_field_mono_nfchoa_sht(xrange,yrange,zrange, Dlsnm, f, conf);
% compute fields from shifted driving functions
Ppw_shift = sound_field_mono_nfchoa_sht(xrange,yrange,zrange, Dpwnm_shift, f, conf);
Pps_shift = sound_field_mono_nfchoa_sht(xrange,yrange,zrange, Dpsnm_shift, f, conf);
Pls_shift = sound_field_mono_nfchoa_sht(xrange,yrange,zrange, Dlsnm_shift, f, conf);
% plot
plot_sound_field(Ppw, x1,y1,z1, [], conf);
title('2.5D NFCHOA (sht domain): plane wave');
plot_sound_field(Pps, x1,y1,z1, [], conf);
title('2.5D NFCHOA (sht domain): point source');
plot_sound_field(Pls, x1,y1,z1, [], conf);
title('2.5D NFCHOA (sht domain): line source');
plot_sound_field(Ppw_shift, x1,y1,z1, [], conf);
title('2.5D NFCHOA (spatial domain): plane wave (shifted expansion)');
plot_sound_field(Pps_shift, x1,y1,z1, [], conf);
title('2.5D NFCHOA (spatial domain): point source (shifted expansion)');
plot_sound_field(Pls_shift, x1,y1,z1, [], conf);
title('2.5D NFCHOA (spatial domain): line source (shifted expansion)');

%% generic 3D NFCHOA in spherical harmonics domain
conf.dimension = '3D';
% compute sht of driving functions
Dpwnm = driving_function_mono_nfchoa_sht_sphexp(Apwnm_3D,f,conf);
Dpsnm = driving_function_mono_nfchoa_sht_sphexp(Apsnm_3D,f,conf);
Dlsnm = driving_function_mono_nfchoa_sht_sphexp(Alsnm_3D,f,conf);
% compute fields
Ppw = sound_field_mono_nfchoa_sht(xrange,yrange,zrange, Dpwnm, f, conf);
Pps = sound_field_mono_nfchoa_sht(xrange,yrange,zrange, Dpsnm, f, conf);
Pls = sound_field_mono_nfchoa_sht(xrange,yrange,zrange, Dlsnm, f, conf);
% plot
plot_sound_field(Ppw, x1,y1,z1, [], conf);
title('3D NFCHOA (sht domain): plane wave');
plot_sound_field(Pps, x1,y1,z1, [], conf);
title('3D NFCHOA (sht domain): point source');
plot_sound_field(Pls, x1,y1,z1, [], conf);
title('3D NFCHOA (sht domain): line source');