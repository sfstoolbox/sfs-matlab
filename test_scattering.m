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
xs = [0,  0.5, 0];  % position of point source
f = 1000;
xrange = [-2 2];
yrange = [-2 2];
zrange = 0;

% scatterer
sigma = inf;  % admittance of scatterer (inf to soft scatterer)
R = 0.3;
xq = [ 0, -1.0, 0];
xt = [ 0.5, 0.5, 0]*R;
conf.xref = xq;
     
display(conf.scattering)

%% Spherical Expansion
% spherical expansion
A1sph = sphexp_mono_pw(ns,f,xq,conf);
A1sph_shift = sphexp_mono_pw(ns,f,xq+xt,conf);
A2sph = sphexp_mono_ps(xs,'R',f,xq,conf);
A2sph_shift = sphexp_mono_ps(xs,'R', f,xq+xt,conf);
% regular-to-regular spherical reexpansion (translatory shift)
RRsph = sphexp_mono_translation(-xt, 'RR', f, conf);
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
% loudspeakers
x0 = secondary_source_positions(conf);

% compute driving functions
D1sph = driving_function_mono_wfs_sphexp(x0(:,1:3),x0(:,4:6),A1sph,'R',f,xq,conf);
D2sph = driving_function_mono_wfs_sphexp(x0(:,1:3),x0(:,4:6),A2sph,'R',f,xq,conf);

% compute fields
P1sphwfs = sound_field_mono(xrange,yrange,zrange,x0,'ps',D1sph,f,conf);
P2sphwfs = sound_field_mono(xrange,yrange,zrange,x0,'ps',D2sph,f,conf);

% plot
[~,~,~,x1,y1,z1] = xyz_grid(xrange,yrange,zrange,conf);

plot_sound_field(P1sphwfs ,x1,y1,z1, x0, conf);
title('plane wave');
plot_sound_field(P2sphwfs ,x1,y1,z1, x0, conf);
title('point source');

%% WFS Reproduction of focused source using time reversal

% singular spherical expansion of point source
B1sph = sphexp_mono_ps(xs, 'S', f, xs, conf);

% compute driving functions
D3sph = driving_function_mono_wfs_sphexp(x0(:,1:3),x0(:,4:6),B1sph,'S',f,xs,conf);

% compute fields (point and focused source)
P3sphwfs = sound_field_mono(xrange,yrange,zrange,x0,'ps',D3sph,f,conf);
P3sphwfs_conj = sound_field_mono(xrange,yrange,zrange,x0,'ps',conj(D3sph),f,conf);

% plot
plot_sound_field(P3sphwfs ,x1,y1,z1, x0, conf);
title('point source');
plot_sound_field(P3sphwfs_conj ,x1,y1,z1, x0, conf);
title('focused source');

%% WFS Reproduction using virtual Scatterer and time reversal
% loudspeakers
x0 = secondary_source_positions(conf);

% compute timereversed incident field
A1sph_timereversed = sphexp_mono_timereverse(A1sph);
A2sph_timereversed = sphexp_mono_timereverse(A2sph);
% compute scattered field
B1sph = sphexp_mono_scatter(A1sph_timereversed, R, sigma, f, conf); 
B2sph = sphexp_mono_scatter(A2sph_timereversed, R, sigma, f, conf); 
% compute driving functions
D1sph = driving_function_mono_wfs_sphexp(x0(:,1:3),x0(:,4:6),B1sph,'S',f,xq,conf);
D2sph = driving_function_mono_wfs_sphexp(x0(:,1:3),x0(:,4:6),B2sph,'S',f,xq,conf);
% compute fields
P1sphwfs = sound_field_mono(xrange,yrange,zrange,x0,'ps',conj(D1sph),f,conf);
P2sphwfs = sound_field_mono(xrange,yrange,zrange,x0,'ps',conj(D2sph),f,conf);

% plot
[~,~,~,x1,y1,z1] = xyz_grid(xrange,yrange,zrange,conf);

plot_sound_field(P1sphwfs ,x1,y1,z1, x0, conf);
title('plane wave');
plot_scatterer(xq,R);
plot_sound_field(P2sphwfs ,x1,y1,z1, x0, conf);
title('point source');
plot_scatterer(xq,R);

%% generic NFCHOA
% loudspeakers
x0 = secondary_source_positions(conf);
% compute driving functions
D1sph = driving_function_mono_nfchoa_sphexp(x0(:,1:3), A1sph, f, conf);
D2sph = driving_function_mono_nfchoa_sphexp(x0(:,1:3), A2sph, f, conf);
% compute fields
P1sphhoa = sound_field_mono(xrange,yrange,zrange,x0,'ps',D1sph,f,conf);
P2sphhoa = sound_field_mono(xrange,yrange,zrange,x0,'ps',D2sph,f,conf);

% plot
[~,~,~,x1,y1,z1] = xyz_grid(xrange,yrange,zrange,conf);

plot_sound_field(P1sphhoa, x1,y1,z1, x0, conf);
title('plane wave');
plot_sound_field(P2sphhoa ,x1,y1,z1, x0, conf);
title('point source');

%% Cylindrical expansion
% A1cyl = cylexpR_mono_pw(ns,f,xq,conf);  % regular expansion plane wave
% A2cyl = cylexpR_mono_pw(ns,f,xq+xt,conf);  % regular expansion plane wave
% % % regular-to-regular cylindrical reexpansion (translatory shift)
%[A1cylt, RR1cyl] = cylexpRR_mono(A1cyl, -xt, f, conf);

% %% Scattering
% % scattering with single sphere
% B1sph = sphexpS_mono_scatter(A1sph, R, sigma, f, conf);  
% B2sph = sphexpS_mono_scatter(A2sph, R, sigma, f, conf);
% % scattering with single cylinder
% B1cyl = cylexpS_mono_scatter(A1cyl, R, sigma, f, conf);
% % singular-to-regular cylindrical reexpansion (translatory shift)
% [B1cylt, SR12cyl] = cylexpSR_mono(B1cyl, -xt, f, conf);
% 
% % multiple scattering
% Bmcyl = cylexpS_mono_multiscatter([A1cyl, A2cyl], [xq; xq + xt], R, ...
%   sigma, f, conf);
% 
% %% WFS Driving Functions
% % driving functions for scattering with single sphere
% x0 = secondary_source_positions(conf);
% x0 = secondary_source_selection(x0,ns,'pw');
% x0 = secondary_source_tapering(x0,conf);
% D1sph = driving_function_mono_wfs_sphexpS(x0(:,1:3),x0(:,4:6),B1sph,f,xq,conf);
% 
% x0 = secondary_source_positions(conf);
% x0 = secondary_source_selection(x0,xs,'ps');
% x0 = secondary_source_tapering(x0,conf);
% D2sph = driving_function_mono_wfs_sphexpS(x0(:,1:3),x0(:,4:6),B2sph,f,xq,conf);
% 
% % driving functions for scattering with single cylinder
% x0 = secondary_source_positions(conf);
% x0 = secondary_source_selection(x0,ns,'pw');
% x0 = secondary_source_tapering(x0,conf);
% D1cyl = driving_function_mono_wfs_cylexpS(x0(:,1:3),x0(:,4:6),B1cyl,f,xq,conf);
% 
% % driving functions for scattering with two cylinders
% x0 = secondary_source_positions(conf);
% x0 = secondary_source_selection(x0,ns,'pw');
% x0 = secondary_source_tapering(x0,conf);
% Dm1cyl = driving_function_mono_wfs_cylexpS(x0(:,1:3),x0(:,4:6),Bmcyl(:,1),f,xq,conf);
% Dm2cyl = driving_function_mono_wfs_cylexpS(x0(:,1:3),x0(:,4:6),Bmcyl(:,2),f,xq + xt,conf);
% 
% %% WFS Sound Fields + Plotting
% % scattering with single sphere
% P1sphwfs = sound_field_mono(xrange,yrange,zrange,x0,'ps',D1sph,f,conf);
% P2sphwfs = sound_field_mono(xrange,yrange,zrange,x0,'ps',D2sph,f,conf);
% % scattering with single cylinder
% P1cylwfs = sound_field_mono(xrange,yrange,zrange,x0,'ls',D1cyl,f,conf);
% % scattering with two cylinders
% Pmcylwfs = sound_field_mono(xrange,yrange,zrange,x0,'ls', Dm1cyl + Dm2cyl,f,conf);
% 
% [~,~,~,x1,y1,z1] = xyz_grid(xrange,yrange,zrange,conf);
% 
% % scattering with single sphere
% plot_sound_field(P1sphwfs ,x1,y1,z1, x0, conf);
% plot_scatterer(xq,R);
% title('WFS: scattering with single sphere (pw)');
% plot_sound_field(P2sphwfs ,x1,y1,z1, x0, conf);
% plot_scatterer(xq,R);
% title('WFS: scattering with single sphere (ps)');
% % scattering with single cylinder
% plot_sound_field(P1cylwfs ,x1,y1,z1, x0, conf);
% plot_scatterer(xq,R);
% title('scattering with single cylinder (pw)');
% % scattering with two cylinders
% plot_sound_field(Pmcylwfs ,x1,y1,z1, x0, conf);
% plot_scatterer(xq,R);
% title('scattering with single cylinder (pw)');
% 

% 
% [Jcyln, Hcyln, Ycyln] = ...
%   cylbasis_mono_XYZgrid(xrange,yrange,zrange,f,xq,conf);
% % 
% [Jcylnt, Hcylnt, Ycylnt] = ...
%   cylbasis_mono_XYZgrid(xrange,yrange,zrange,f,xq+xt,conf);
%% Spherical Expansion Sound Fields + Plotting
% incident fields 

% P1cyl = sound_field_mono_basis(A1cyl, Jcyln, Ycyln, conf);
% P1cylt = sound_field_mono_basis(A1cylt, Jcylnt, Ycylnt, conf);
% P2cyl = sound_field_mono_basis(A2cyl, Jcylnt, Ycylnt, conf);

% % scattering with single sphere
% P1sphscat = sound_field_mono_basis(B1sph, Hsphn, Ysphnm, conf);
% P2sphscat = sound_field_mono_basis(B2sph, Hsphn, Ysphnm, conf);
% % scattering with single cylinder
% P1cylscat = sound_field_mono_basis(B1cyl, Hcyln, Ycyln, conf);
% P2cylscat = sound_field_mono_basis(B1cylt, Jcylnt, Ycylnt, conf);
% % multiple scattering with two cylinders
% P1cylmscat = sound_field_mono_basis(Bmcyl(:,1), Hcyln, Ycyln, conf);
% P2cylmscat = sound_field_mono_basis(Bmcyl(:,2), Hcylnt, Ycylnt, conf);

%[~,~,~,x1,y1,z1] = xyz_grid(xrange,yrange,zrange,conf);

% scattering with single sphere

% plot_sound_field(P1sphscat ,x1,y1,z1, [], conf);
% plot_scatterer(xq,R);
% title('scattered field');
% plot_sound_field(P1sph + P1sphscat ,x1,y1,z1, [], conf);
% plot_scatterer(xq,R);
% title('incident field + scattered field');

% scattering with single sphere
% plot_sound_field(P2sph ,x1,y1,z1, [], conf);
% plot_scatterer(xq,R);
% title('incident field');
% plot_sound_field(P2sphscat ,x1,y1,z1, [], conf);
% plot_scatterer(xq,R);
% title('scattered field');
% plot_sound_field(P2sph + P2sphscat ,x1,y1,z1, [], conf);
% plot_scatterer(xq,R);
% title('incident field + scattered field');

% % scattering with single cylinder
% plot_sound_field(P1cyl ,x1,y1,z1, [], conf);
% plot_scatterer(xq,R);
% title('incident field 1');
% plot_sound_field(P1cylt ,x1,y1,z1, [], conf);
% plot_scatterer(xq + xt,R);
% title('incident field 1 (shifted)');
% plot_sound_field(P2cyl ,x1,y1,z1, [], conf);
% plot_scatterer(xq + xt,R);
% title('incident field 2');
% plot_sound_field(P1cylscat ,x1,y1,z1, [], conf);
% plot_scatterer(xq,R);
% title('scattered field 1');
% plot_sound_field(P2cylscat ,x1,y1,z1, [], conf);
% plot_scatterer(xq + xt,R);
% title('scattered field (shifted) 1');
% plot_sound_field(P1cyl + P1cylscat ,x1,y1,z1, [], conf);
% plot_scatterer(xq,R);
% title('incident field + scattered field 1');
% 
% % scattering with two cylinders
% plot_sound_field(P1cylmscat ,x1,y1,z1, [], conf);
% plot_scatterer(xq,R);
% title('two cylinders: scattered field 1');
% plot_sound_field(P1cyl + P1cylmscat ,x1,y1,z1, [], conf);
% plot_scatterer(xq,R);
% title('two cylinders: scattered field + incident field 1');
% plot_sound_field(P2cylmscat,x1,y1,z1, [], conf);
% plot_scatterer(xq + xt,R);
% title('two cylinders: scattered field 2');
% plot_sound_field(P2cyl + P2cylmscat ,x1,y1,z1, [], conf);
% plot_scatterer(xq + xt,R);
% title('two cylinders: scattered field + incident field 2');
% plot_sound_field(P1cyl + P1cylmscat + P2cyl + P2cylmscat ,x1,y1,z1, [], conf);
% plot_scatterer(xq + xt,R);
% title('two cylinders: scattered field + incident field all');