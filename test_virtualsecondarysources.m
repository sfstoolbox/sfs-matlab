%% initialize
close all;
clear variables;
% SFS Toolbox
addpath('~/projects/sfstoolbox'); SFS_start;

%% Parameters
conf = SFS_config_example;
conf.plot.useplot = false;
conf.showprogress = true;
conf.resolution = 400;
conf.plot.loudspeakers = true;
conf.plot.realloudspeakers = true;

% config for virtual array
conf.virtual_secondary_sources.size = 0.4;
conf.virtual_secondary_sources.center = [0, 0.5, 0];
conf.virtual_secondary_sources.geometry = 'circular';
conf.virtual_secondary_sources.number = 20;

% config for real array
conf.dimension = '3D';
conf.secondary_sources.geometry = 'circular';
conf.secondary_sources.number = 56;
conf.secondary_sources.size = 3;
conf.secondary_sources.center = [0, 0, 0];
conf.xref = conf.virtual_secondary_sources.center;

xs = [0, -1, 0];  % propagation direction of plane wave
src = 'pw'; 
f = 4000;
xrange = [-2 2];
yrange = [-2 2];
zrange = 0;

%%
x0 = secondary_source_positions(conf);

[D, xv, x0] = driving_function_mono_localwfs(x0,xs,src,f,conf);

[P, x1, y1, z1] = sound_field_mono(xrange,yrange,zrange,x0,'ps',D,f,conf);

dimensions = xyz_axes_selection(x1,y1,z1);

plot_sound_field(P,x1,y1,z1, x0, conf);
hold on
  draw_loudspeakers(xv,dimensions,conf);
hold off
