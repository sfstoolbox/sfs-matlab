%% initialize
% close all;
clear variables;
% SFS Toolbox
SFS_start;

%% Parameters
conf = SFS_config;
conf.showprogress = true;

% plotting
conf.plot.usedb = false;
conf.plot.useplot = false;
conf.usenormalisation = true;
conf.resolution = 200;

X = [-2 2];
Y = [-2 2];
Z = 0;

f = 1000;
ns = [0, -1, 0];  % propagation direction of plane wave
xs = [0,  3, 0];  % position of point source

xt = [ 0.5, 0.5, 0];  % position of the sweet spot
rt = 0.3;  % "size" of the sweet spot
Ncet = circexp_truncation_order(norm(xt)+rt, f, 1e-6, conf);  % 2.5DHOA order for translation
Nce = circexp_truncation_order(rt, f, 1e-6, conf);  % 2.5DHOA order for sweet spot
Nse = sphexp_truncation_order(rt, f, 1e-6, conf);  % 3DHOA order for sweet spot

% secondary sources
conf.secondary_sources.geometry = 'circular';
conf.secondary_sources.number = Nce*2+1;
conf.secondary_sources.size = 3;
conf.secondary_sources.center = [0, 0, 0];

xq = conf.secondary_sources.center;  % expansion center
conf.xref = xq;  % reference position

%% Spherical Basis Functions
[jn, ~, Ynm] = sphbasis_mono_grid(X,Y,Z,max(Ncet,Nse),f,xq,conf);

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
% computed truncated expansion
Apwnm_shift_trunc = sphexp_truncate(Apwnm_shift, Nce, -round(2*pi*f/conf.c*norm(xt)*sind(90-atan2d(-xt(2),-xt(1)))));
Apsnm_shift_trunc = sphexp_truncate(Apsnm_shift, Nce, -round(2*pi*f/conf.c*norm(xt)*sind(90-atan2d(-xt(2),-xt(1)))));
Alsnm_shift_trunc = sphexp_truncate(Alsnm_shift, Nce, -round(2*pi*f/conf.c*norm(xt)*sind(90-atan2d(-xt(2),-xt(1)))));

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
% compute driving functions from shifted fields
Dpw_shift_trunc = driving_function_mono_nfchoa_sphexp(x0(:,1:3), Apwnm_shift_trunc, f, conf);
Dps_shift_trunc = driving_function_mono_nfchoa_sphexp(x0(:,1:3), Apsnm_shift_trunc, f, conf);
Dls_shift_trunc = driving_function_mono_nfchoa_sphexp(x0(:,1:3), Alsnm_shift_trunc, f, conf);
% compute fields
Ppw = sound_field_mono(X,Y,Z,x0,'ps',Dpw,f,conf);
Pps = sound_field_mono(X,Y,Z,x0,'ps',Dps,f,conf);
Pls = sound_field_mono(X,Y,Z,x0,'ls',Dls,f,conf);
% compute fields from shifted driving functions
Ppw_shift = sound_field_mono(X,Y,Z,x0,'ps',Dpw_shift,f,conf);
Pps_shift = sound_field_mono(X,Y,Z,x0,'ps',Dps_shift,f,conf);
Pls_shift = sound_field_mono(X,Y,Z,x0,'ls',Dls_shift,f,conf);
% compute fields from shifted + bandlimited driving functions
Ppw_shift_trunc = sound_field_mono(X,Y,Z,x0,'ps',Dpw_shift_trunc,f,conf);
Pps_shift_trunc = sound_field_mono(X,Y,Z,x0,'ps',Dps_shift_trunc,f,conf);
Pls_shift_trunc = sound_field_mono(X,Y,Z,x0,'ls',Dls_shift_trunc,f,conf);
% plot

plot_sound_field(Ppw, X, Y, Z, [], conf);
title('2.5D NFCHOA (spatial domain): plane wave');
plot_sound_field(Pps, X, Y, Z, [], conf);
title('2.5D NFCHOA (spatial domain): point source');
plot_sound_field(Pls, X, Y, Z, [], conf);
title('2.5D NFCHOA (spatial domain): line source');
plot_sound_field(Ppw_shift, X, Y, Z, [], conf);
title('2.5D NFCHOA (spatial domain): plane wave (shifted expansion)');
plot_sound_field(Pps_shift, X, Y, Z, [], conf);
title('2.5D NFCHOA (spatial domain): point source (shifted expansion)');
plot_sound_field(Pls_shift, X, Y, Z, [], conf);
title('2.5D NFCHOA (spatial domain): line source (shifted expansion)');
plot_sound_field(Ppw_shift_trunc, X, Y, Z, [], conf);
title('2.5D NFCHOA (spatial domain): plane wave (shifted + truncated expansion)');
plot_sound_field(Pps_shift_trunc, X, Y, Z, [], conf);
title('2.5D NFCHOA (spatial domain): point source (shifted + truncated expansion)');
plot_sound_field(Pls_shift_trunc, X, Y, Z, [], conf);
title('2.5D NFCHOA (spatial domain): line source (shifted + truncated expansion)');

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
% compute sht of driving functions from shifted fields
Dpwnm_shift_trunc = driving_function_mono_nfchoa_sht_sphexp(Apwnm_shift_trunc, f, conf);
Dpsnm_shift_trunc = driving_function_mono_nfchoa_sht_sphexp(Apsnm_shift_trunc, f, conf);
Dlsnm_shift_trunc = driving_function_mono_nfchoa_sht_sphexp(Alsnm_shift_trunc, f, conf);
% compute spherical expansion of reproduced sound field
Ppwnm = sphexp_mono_sht(Dpwnm,'R',f,conf);
Ppsnm = sphexp_mono_sht(Dpsnm,'R',f,conf);
Plsnm = sphexp_mono_sht(Dlsnm,'R',f,conf);
% compute spherical expansion of reproduced sound field from shifted driving functions
Ppwnm_shift = sphexp_mono_sht(Dpwnm_shift,'R',f,conf);
Ppsnm_shift = sphexp_mono_sht(Dpsnm_shift,'R',f,conf);
Plsnm_shift = sphexp_mono_sht(Dlsnm_shift,'R',f,conf);
% compute spherical expansion of reproduced sound field from shifted and truncated driving functions
Ppwnm_shift_trunc = sphexp_mono_sht(Dpwnm_shift_trunc,'R',f,conf);
Ppsnm_shift_trunc = sphexp_mono_sht(Dpsnm_shift_trunc,'R',f,conf);
Plsnm_shift_trunc = sphexp_mono_sht(Dlsnm_shift_trunc,'R',f,conf);
% compute fields
Ppw = sound_field_mono_sphbasis(Ppwnm, jn, Ynm);
Pps = sound_field_mono_sphbasis(Ppsnm, jn, Ynm);
Pls = sound_field_mono_sphbasis(Plsnm, jn, Ynm);
% compute fields from shifted driving functions
Ppw_shift = sound_field_mono_sphbasis(Ppwnm_shift, jn, Ynm);
Pps_shift = sound_field_mono_sphbasis(Ppsnm_shift, jn, Ynm);
Pls_shift = sound_field_mono_sphbasis(Plsnm_shift, jn, Ynm);
% compute fields from shifted driving functions
Ppw_shift_trunc = sound_field_mono_sphbasis(Ppwnm_shift_trunc, jn, Ynm);
Pps_shift_trunc = sound_field_mono_sphbasis(Ppsnm_shift_trunc, jn, Ynm);
Pls_shift_trunc = sound_field_mono_sphbasis(Plsnm_shift_trunc, jn, Ynm);
% plot

plot_sound_field(Ppw, X, Y, Z, [], conf);
title('2.5D NFCHOA (sht domain): plane wave');
plot_sound_field(Pps, X, Y, Z, [], conf);
title('2.5D NFCHOA (sht domain): point source');
plot_sound_field(Pls, X, Y, Z, [], conf);
title('2.5D NFCHOA (sht domain): line source');
plot_sound_field(Ppw_shift, X, Y, Z, [], conf);
title('2.5D NFCHOA (spatial domain): plane wave (shifted expansion)');
plot_sound_field(Pps_shift, X, Y, Z, [], conf);
title('2.5D NFCHOA (spatial domain): point source (shifted expansion)');
plot_sound_field(Pls_shift, X, Y, Z, [], conf);
title('2.5D NFCHOA (spatial domain): line source (shifted expansion)');
plot_sound_field(Ppw_shift_trunc, X, Y, Z, [], conf);
title('2.5D NFCHOA (spatial domain): plane wave (shifted + truncated expansion)');
plot_sound_field(Pps_shift_trunc, X, Y, Z, [], conf);
title('2.5D NFCHOA (spatial domain): point source (shifted + truncated expansion)');
plot_sound_field(Pls_shift_trunc, X, Y, Z, [], conf);
title('2.5D NFCHOA (spatial domain): line source (shifted + truncated expansion)');

%% generic 3D NFCHOA in spherical harmonics domain
conf.dimension = '3D';
% compute sht of driving functions
Dpwnm = driving_function_mono_nfchoa_sht_sphexp(Apwnm_3D,f,conf);
Dpsnm = driving_function_mono_nfchoa_sht_sphexp(Apsnm_3D,f,conf);
Dlsnm = driving_function_mono_nfchoa_sht_sphexp(Alsnm_3D,f,conf);
% compute spherical expansion of reproduced sound field
Ppwnm = sphexp_mono_sht(Dpwnm,'R',f,conf);
Ppsnm = sphexp_mono_sht(Dpsnm,'R',f,conf);
Plsnm = sphexp_mono_sht(Dlsnm,'R',f,conf);
% compute fields
Ppw = sound_field_mono_sphbasis(Ppwnm, jn, Ynm);
Pps = sound_field_mono_sphbasis(Ppsnm, jn, Ynm);
Pls = sound_field_mono_sphbasis(Plsnm, jn, Ynm);
% plot

plot_sound_field(Ppw, X, Y, Z, [], conf);
title('3D NFCHOA (sht domain): plane wave');
plot_sound_field(Pps, X, Y, Z, [], conf);
title('3D NFCHOA (sht domain): point source');
plot_sound_field(Pls, X, Y, Z, [], conf);
title('3D NFCHOA (sht domain): line source');