%% initialize
close all;
clear variables;
% SFS Toolbox
addpath('~/projects/sfstoolbox'); SFS_start;

%% Parameters
conf = SFS_config_example;

conf.dimension = '3D';
conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.number = 40;
conf.secondary_sources.size = 6;
conf.secondary_sources.center = [0, 2, 0];

conf.plot.useplot = false;

conf.scattering.Nse = 23;
conf.scattering.Nce = 23;
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
xq = [ 2*R, -R, 0];
xt = [ -3*R, 0, 0];
conf.xref = xq;
     
%% Expansion
% spherical expansion
A1sph = sphexpR_mono_pw(ns,f,xq,conf);  % regular expansion plane wave
A2sph = sphexpR_mono_ps(xs,f,xq,conf);  % regular expansion point source
% cylindrical expansion
A1cyl = cylexpR_mono_pw(ns,f,xq,conf);  % regular expansion plane wave
A2cyl = cylexpR_mono_pw(ns,f,xq+xt,conf);  % regular expansion plane wave
% regular-to-regular cylindrical reexpansion (translatory shift)
[A1cylt, RR1cyl] = cylexpRR_mono(A1cyl, -xt, f, conf);

%% Scattering
% scattering with single sphere
B1sph = sphexpS_mono_scatter(A1sph, R, sigma, f, conf);  
B2sph = sphexpS_mono_scatter(A2sph, R, sigma, f, conf);
% scattering with single cylinder
B1cyl = cylexpS_mono_scatter(A1cyl, R, sigma, f, conf);
% singular-to-regular cylindrical reexpansion (translatory shift)
[B1cylt, SR12cyl] = cylexpSR_mono(B1cyl, -xt, f, conf);

% multiple scattering
Bmcyl = cylexpS_mono_multiscatter([A1cyl, A2cyl], [xq; xq + xt], R, ...
  sigma, f, conf);

%% WFS Driving Functions
% driving functions for scattering with single sphere
x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,ns,'pw');
x0 = secondary_source_tapering(x0,conf);
D1sph = driving_function_mono_wfs_sphexpS(x0(:,1:3),x0(:,4:6),B1sph,f,xq,conf);

x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,xs,'ps');
x0 = secondary_source_tapering(x0,conf);
D2sph = driving_function_mono_wfs_sphexpS(x0(:,1:3),x0(:,4:6),B2sph,f,xq,conf);

% driving functions for scattering with single cylinder
x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,ns,'pw');
x0 = secondary_source_tapering(x0,conf);
D1cyl = driving_function_mono_wfs_cylexpS(x0(:,1:3),x0(:,4:6),B1cyl,f,xq,conf);

% driving functions for scattering with two cylinders
x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,ns,'pw');
x0 = secondary_source_tapering(x0,conf);
Dm1cyl = driving_function_mono_wfs_cylexpS(x0(:,1:3),x0(:,4:6),Bmcyl(:,1),f,xq,conf);
Dm2cyl = driving_function_mono_wfs_cylexpS(x0(:,1:3),x0(:,4:6),Bmcyl(:,2),f,xq + xt,conf);

%% WFS Sound Fields + Plotting
% scattering with single sphere
P1sphwfs = sound_field_mono(xrange,yrange,zrange,x0,'ps',D1sph,f,conf);
P2sphwfs = sound_field_mono(xrange,yrange,zrange,x0,'ps',D2sph,f,conf);
% scattering with single cylinder
P1cylwfs = sound_field_mono(xrange,yrange,zrange,x0,'ls',D1cyl,f,conf);
% scattering with two cylinders
Pmcylwfs = sound_field_mono(xrange,yrange,zrange,x0,'ls', Dm1cyl + Dm2cyl,f,conf);

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
% scattering with two cylinders
plot_sound_field(Pmcylwfs ,x1,y1,z1, x0, conf);
plot_scatterer(xq,R);
title('scattering with single cylinder (pw)');

%% Evaluate spherical and cylindrical basis functions 
[Jsphn, Hsphn, Ysphnm] = ...
  sphbasis_mono_XYZgrid(xrange,yrange,zrange,f,xq,conf);

[Jcyln, Hcyln, Ycyln] = ...
  cylbasis_mono_XYZgrid(xrange,yrange,zrange,f,xq,conf);

[Jcylnt, Hcylnt, Ycylnt] = ...
  cylbasis_mono_XYZgrid(xrange,yrange,zrange,f,xq+xt,conf);
%% Spherical Expansion Sound Fields + Plotting
% incident fields 
P1sph = sound_field_mono_basis(A1sph, Jsphn, Ysphnm, conf);
P2sph = sound_field_mono_basis(A2sph, Jsphn, Ysphnm, conf);
P1cyl = sound_field_mono_basis(A1cyl, Jcyln, Ycyln, conf);
P1cylt = sound_field_mono_basis(A1cylt, Jcylnt, Ycylnt, conf);
P2cyl = sound_field_mono_basis(A2cyl, Jcylnt, Ycylnt, conf);

% scattering with single sphere
P1sphscat = sound_field_mono_basis(B1sph, Hsphn, Ysphnm, conf);
P2sphscat = sound_field_mono_basis(B2sph, Hsphn, Ysphnm, conf);
% scattering with single cylinder
P1cylscat = sound_field_mono_basis(B1cyl, Hcyln, Ycyln, conf);
P2cylscat = sound_field_mono_basis(B1cylt, Jcylnt, Ycylnt, conf);
% multiple scattering with two cylinders
P1cylmscat = sound_field_mono_basis(Bmcyl(:,1), Hcyln, Ycyln, conf);
P2cylmscat = sound_field_mono_basis(Bmcyl(:,2), Hcylnt, Ycylnt, conf);

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
title('incident field 1');
plot_sound_field(P1cylt ,x1,y1,z1, [], conf);
plot_scatterer(xq + xt,R);
title('incident field 1 (shifted)');
plot_sound_field(P2cyl ,x1,y1,z1, [], conf);
plot_scatterer(xq + xt,R);
title('incident field 2');
plot_sound_field(P1cylscat ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('scattered field 1');
plot_sound_field(P2cylscat ,x1,y1,z1, [], conf);
plot_scatterer(xq + xt,R);
title('scattered field (shifted) 1');
plot_sound_field(P1cyl + P1cylscat ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('incident field + scattered field 1');

% scattering with two cylinders
plot_sound_field(P1cylmscat ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('two cylinders: scattered field 1');
plot_sound_field(P1cyl + P1cylmscat ,x1,y1,z1, [], conf);
plot_scatterer(xq,R);
title('two cylinders: scattered field + incident field 1');
plot_sound_field(P2cylmscat,x1,y1,z1, [], conf);
plot_scatterer(xq + xt,R);
title('two cylinders: scattered field 2');
plot_sound_field(P2cyl + P2cylmscat ,x1,y1,z1, [], conf);
plot_scatterer(xq + xt,R);
title('two cylinders: scattered field + incident field 2');
plot_sound_field(P1cyl + P1cylmscat + P2cyl + P2cylmscat ,x1,y1,z1, [], conf);
plot_scatterer(xq + xt,R);
title('two cylinders: scattered field + incident field all');