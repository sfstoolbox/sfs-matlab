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
conf.plot.usedb = true;
conf.tapwinlen = 1.0;

% config for virtual array
conf.localsfs.size = 0.4;
conf.localsfs.center = [0.0, 0, 0];

conf.localsfs.vss.geometry = 'circular';
conf.localsfs.vss.number = 500;
conf.localsfs.vss.sampling = 'equi';
conf.localsfs.vss.logratio = 1.0;
conf.localsfs.vss.method = 'wfs';
conf.localsfs.vss.type = 'ls';
conf.localsfs.vss.consider_local_area = true;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
conf.localsfs.vss.usetapwin = true;
conf.localsfs.vss.tapwinlen = 0.3;

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
xs = [0.0, -1.0, 0];  % propagation direction of plane wave
src = 'pw';
xrange = [-conf.secondary_sources.size/2, conf.secondary_sources.size/2];
yrange = [-conf.secondary_sources.size/2, conf.secondary_sources.size/2];
zrange = 0;

x0 = secondary_source_positions(conf);
x0s = secondary_source_selection(x0, xs, src);
x0s = secondary_source_tapering(x0s, conf);

[d1, xactive, xv] = driving_function_imp_localwfs(x0s,xs,src,conf);

d2 = driving_function_imp_wfs(x0s,xs,src,conf);

%%
[p,x,y,z] = sound_field_imp(xrange,yrange,zrange, x0s, 'ls', d2, 190,conf);
plot_sound_field(p,x,y,z,x0s,conf);
[p,x,y,z] = sound_field_imp(xrange,yrange,zrange,xactive, 'ls',d1, 150,conf);
plot_sound_field(p,x,y,z,xactive,conf);
hold on
  draw_loudspeakers(xv, [1 1 0], conf);
hold off
