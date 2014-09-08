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
conf.localsfs.size = 0.6;
conf.localsfs.center = [0, 0, 0];
conf.localsfs.geometry = 'linear';

conf.localsfs.vss.number = 16;
conf.localsfs.vss.sampling = 'equi';
conf.localsfs.vss.logratio = 1.0;
conf.localsfs.vss.method = 'wfs';
conf.localsfs.vss.type = 'ls';
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
conf.localsfs.vss.tapering = false;

% config for real array
conf.dimension = '2D';
conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.number = 16;
conf.secondary_sources.size = 2;
conf.secondary_sources.center = [0, 2, 0];
conf.driving_functions = 'default';
conf.xref = conf.localsfs.center;

X = conf.localsfs.center;
xs = [0.0, -1.0, 0];  % propagation direction of plane wave
src = 'pw';
xrange = [-conf.localsfs.size/2, conf.localsfs.size/2] + conf.xref(1);
yrange = [-conf.localsfs.size/2, conf.localsfs.size/2] + conf.xref(2);
zrange = 0;
f = 9000;

x0 = secondary_source_positions(conf);

[D, xactive, xv] = driving_function_mono_localwfs(x0,xs,src,f,conf);
[P, x, y, z] = sound_field_mono(xrange,yrange,zrange,xactive,'ls',D,f,conf);
plot_sound_field(P,x,y,z,x0,conf);
hold on
draw_loudspeakers(xv,[1 1 0],conf);
hold off
