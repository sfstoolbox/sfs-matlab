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
conf.localsfs.size = 0.4;
conf.localsfs.center = [0.5, 0, 0];
conf.localsfs.geometry = 'circular';
conf.localsfs.number = 56;
conf.localsfs.method = 'wfs';

% config for real array
conf.dimension = '2.5D';
conf.secondary_sources.geometry = 'circular';
conf.secondary_sources.number = 56;
conf.secondary_sources.size = 3;
conf.secondary_sources.center = [0, 0, 0];
conf.xref = conf.localsfs.center;

xs = [0.0, 1.0, 0];  % propagation direction of plane wave
src = 'ps'; 
f = 4000;
xrange = [-2 2];
yrange = [-2 2];
zrange = 0;
%%
x0 = secondary_source_positions(conf);

[D, xv, x0] = driving_function_mono_localwfs(x0,xs,src,f,conf);

[P, x1, y1, z1] = sound_field_mono(xrange,yrange,zrange,x0,'ps',D,f,conf);
Pgt = sound_field_mono(xrange,yrange,zrange,[xs, 1 0 0, 1],src,1,f,conf);

dimensions = xyz_axes_selection(x1,y1,z1);

plot_sound_field(P,x1,y1,z1, x0, conf);
hold on
  conf.plot.realloudspeakers = false;
  draw_loudspeakers(xv,dimensions,conf);
hold off
plot_sound_field(Pgt,x1,y1,z1, xv, conf);
conf.plot.usedb = true;
plot_sound_field(abs((P-Pgt)./Pgt),x1,y1,z1, xv, conf);