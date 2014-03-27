%% initialize
close all;
clear variables;
% SFS Toolbox
addpath('~/projects/sfstoolbox'); SFS_start;

%% Parameters
conf = SFS_config_example;

conf.dimension = '2.5D';
conf.driving_functions = 'point_source';
conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.number = 40;
conf.secondary_sources.size = 6;
conf.secondary_sources.center = [0, 2, 0];

conf.plot.useplot = false;

conf.scattering.Nse = 23;
conf.scattering.Nce = 100;
conf.scattering.timereverse = true;  % time reverse wavefield of sph-/cylexpR_

conf.showprogress = true;
conf.resolution = 400;

ns = [0, -1, 0];  % propagation direction of plane wave
xs = [0,  3, 0];  % position of point source
f = 4000;
xrange = [-2 2];
yrange = [-2 2];
zrange = 0;

% scatterer
sigma = inf;  % admittance of scatterer (inf to soft scatterer)
R = 0.3;
xq = [2*R, -R, 0];
conf.xref = xq;
     
%% Expansion
% spherical expansion
A1sph = sphexpR_mono_pw(ns,f,xq,conf);  % regular expansion plane wave
A2sph = sphexpR_mono_ps(xs,f,xq,conf);  % regular expansion point source
% cylindrical expansion
A1cyl = cylexpR_mono_pw(ns,f,xq,conf);  % regular expansion plane wave

%% Scattering
% scattering with single sphere
B1sph = sphexpS_mono_scatter(A1sph, R, sigma, f, conf);  
B2sph = sphexpS_mono_scatter(A2sph, R, sigma, f, conf);
% scattering with single cylinder
B1cyl = cylexpS_mono_scatter(A1cyl, R, sigma, f, conf);

%% WFS Driving Signals
% driving for scattering with single sphere
conf.dimension = '2.5D';
conf.driving_functions = 'point_source';
x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,ns,'pw');
x0 = secondary_source_tapering(x0,conf);
D1sph = driving_function_mono_wfs_sphscatter(x0(:,1:3),x0(:,4:6),B1sph,f,xq,conf);

x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,xs,'ps');
x0 = secondary_source_tapering(x0,conf);
D2sph = driving_function_mono_wfs_sphscatter(x0(:,1:3),x0(:,4:6),B2sph,f,xq,conf);

% driving for scattering with single cylinder
conf.dimension = '3D';
conf.driving_functions = 'line_source';
x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,ns,'pw');
x0 = secondary_source_tapering(x0,conf);
D1cyl = driving_function_mono_wfs_cylscatter(x0(:,1:3),x0(:,4:6),B1cyl,f,xq,conf);

%% WFS Sound Fields + Plotting
% scattering with single sphere
P1sphwfs = sound_field_mono(xrange,yrange,zrange,x0,'ps',D1sph,f,conf);
P2sphwfs = sound_field_mono(xrange,yrange,zrange,x0,'ps',D2sph,f,conf);
% scattering with single cylinder
P1cylwfs = sound_field_mono(xrange,yrange,zrange,x0,'ls',D1cyl,f,conf);

[~,~,~,x1,y1,z1] = xyz_grid(xrange,yrange,zrange,conf);

% scattering with single sphere
plot_sound_field(P1sphwfs ,x1,y1,z1, x0, conf);
plot_scatterer(xq,R);
title('WFS: scattering with single sphere (pw)');
plot_sound_field(P2sphwfs ,x1,y1,z1, x0, conf);
plot_scatterer(xq,R);
title('WFS: scattering with single sphere (ps)');
% scattering with single cylinder
plot_sound_field(P1cylwfs ,x1,y1,z1, x0, conf);
plot_scatterer(xq,R);
title('scattering with single cylinder (pw)');

%% Evaluate spherical and cylindrical basis functions 
[Jsphn, Hsphn, Ysphnm] = ...
  eval_sphbasis_mono_XYZgrid(xrange,yrange,zrange,f,xq,conf);

[Jcyln, Hcyln, Ycyln] = ...
  eval_cylbasis_mono_XYZgrid(xrange,yrange,zrange,f,xq,conf);
%% Spherical Expansion Sound Fields + Plotting
% incident fields 
P1sph = sound_field_mono_sphbasis(A1sph, Jsphn, Ysphnm, conf);
P2sph = sound_field_mono_sphbasis(A2sph, Jsphn, Ysphnm, conf);
P1cyl = sound_field_mono_cylbasis(A1cyl, Jcyln, Ycyln, conf);

% scattering with single sphere
P1sphscat = sound_field_mono_sphbasis(B1sph, Hsphn, Ysphnm, conf);
P2sphscat = sound_field_mono_sphbasis(B2sph, Hsphn, Ysphnm, conf);
% scattering with single cylinder
P1cylscat = sound_field_mono_cylbasis(B1cyl, Hcyln, Ycyln, conf);

[~,~,~,x1,y1,z1] = xyz_grid(xrange,yrange,zrange,conf);

% scattering with single sphere
plot_sound_field(P1sph ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('incident field');
plot_sound_field(P1sphscat ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('scattered field');
plot_sound_field(P1sph + P1sphscat ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('incident field + scattered field');

% scattering with single sphere
plot_sound_field(P2sph ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('incident field');
plot_sound_field(P2sphscat ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('scattered field');
plot_sound_field(P2sph + P2sphscat ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('incident field + scattered field');

% scattering with single cylinder
plot_sound_field(P1cyl ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('incident field');
plot_sound_field(P1cylscat ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('scattered field');
plot_sound_field(P1cyl + P1cylscat ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('incident field + scattered field');


%% Plotting
