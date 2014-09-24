%% initialize
close all;
clear variables;

% startup;

% SFS Toolbox
addpath('~/projects/sfstoolbox'); SFS_start;

%% Parameters
conf = SFS_config_example;
conf.plot.useplot = false;
conf.showprogress = true;
conf.resolution = 400;
conf.plot.loudspeakers = true;
conf.plot.realloudspeakers = false;
conf.usetapwin = true;

% config for virtual array
conf.localsfs.size = 0.4;
conf.localsfs.center = [0.0, 0, 0];

conf.localsfs.vss.geometry = 'circular';
conf.localsfs.vss.number = 14;
conf.localsfs.vss.sampling = 'equi';
conf.localsfs.vss.logratio = 1.0;
conf.localsfs.vss.method = 'wfs';
conf.localsfs.vss.type = 'ls';
conf.localsfs.vss.consider_local_area = true;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
conf.localsfs.vss.tapering = false;

% config for real array
conf.dimension = '2D';
conf.secondary_sources.geometry = 'circular';
conf.secondary_sources.number = 56;
conf.secondary_sources.size = 2;
conf.secondary_sources.center = [0, 0, 0];
conf.driving_functions = 'default';
conf.xref = conf.secondary_sources.center;

conf.wfs.hprefhigh = 6000;

X = conf.localsfs.center;
xs = [0.0, 1.5, 0];  % propagation direction of plane wave
src = 'ps';
xrange = [-1, 1];
yrange = [-1, 1];
zrange = 0;

x0 = secondary_source_positions(conf);

[d1, xactive, xv] = driving_function_imp_localwfs(x0,xs,src,conf);

d2 = driving_function_imp_wfs(x0,xs,src,conf);

subplot(1,2,1), plot(db(easyfft(d1(:,1), conf)));
subplot(1,2,2), plot(db(easyfft(d2(:,1), conf)));
%%
[p,x,y,z] = sound_field_imp(xrange,yrange,zrange,x0, 'ls', d2, 200,conf);
plot_sound_field(p,x,y,z,x0,conf);
[p,x,y,z] = sound_field_imp(xrange,yrange,zrange,xactive, 'ls',d1, 300,conf);
plot_sound_field(p,x,y,z,x0,conf);
hold on
  draw_loudspeakers(xv, [1 1 0], conf);
hold off
