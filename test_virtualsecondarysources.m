%% initialize
close all;
clear variables;

startup;

% SFS Toolbox
addpath('~/projects/sfstoolbox'); SFS_start;

%% Parameters
conf = SFS_config_example;
conf.plot.useplot = false;
conf.showprogress = true;
conf.resolution = 400;
conf.plot.loudspeakers = true;
conf.plot.realloudspeakers = true;
conf.usetapwin = true;

% config for virtual array
conf.localsfs.size = 0.5;
conf.localsfs.center = [0, 0, 0];
conf.localsfs.geometry = 'circular';

conf.localsfs.vss.number = 53;
conf.localsfs.vss.sampling = 'log';
conf.localsfs.vss.method = 'wfs';
conf.localsfs.vss.type = 'ls';
conf.localsfs.vss.ignoress = false;

% config for real array
conf.dimension = '2D';
conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.number = 56;
conf.secondary_sources.size = 6;
conf.secondary_sources.center = [0, 2, 0];
conf.driving_functions = 'default';
conf.xref = conf.localsfs.center;

X = conf.localsfs.center;
xs = [0, -1.0, 0];  % propagation direction of plane wave
src = 'pw'; 
xrange = [-conf.localsfs.size, conf.localsfs.size];
yrange = [-conf.localsfs.size, conf.localsfs.size];
zrange = 0;
f = 9000;

x0 = secondary_source_positions(conf);

%%
for r=[1/10, 1/3, 0.5, 2/3, 1.0, 1.5, 2.0, 3.0, 10]
  conf.localsfs.vss.logratio = r;
  [D, xactive, xv] = driving_function_mono_localwfs(x0,xs,src,f,conf);
  sound_field_mono(xrange,yrange,zrange,xactive,'ls',D,f,conf);
  title(['ratio between shortest and longest distance: ', num2str(r)]);  
  hold on;
  conf.plot.realloudspeakers = false;
  draw_loudspeakers(xv,[1 1 0],conf);
  conf.plot.realloudspeakers = true;
  hold off;
  axis([-conf.localsfs.size, conf.localsfs.size, -conf.localsfs.size, conf.localsfs.size]);
end

%%
for r=[0.5, 1.0, 2]
  conf.localsfs.vss.logratio = r;
  freq_response_localwfs(X,xs,src,conf);
  title(['ratio between shortest and longest distance: ', num2str(r)]);
end