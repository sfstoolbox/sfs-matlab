%% initialize
close all;
clear variables;

%% Parameters
conf = SFS_config_example;
conf.plot.useplot = false;
conf.showprogress = true;
conf.resolution = 400;
conf.plot.loudspeakers = true;
conf.plot.realloudspeakers = false;
conf.usetapwin = true;

% config for virtual array
conf.localsfs.vss.size = 9;
conf.localsfs.vss.center = [0, 0, 0];
conf.localsfs.vss.geometry = 'linear';
conf.localsfs.vss.number = 20;
conf.localsfs.vss.sampling = 'log';
conf.localsfs.method = 'wfs';
% TODO: what does this?
conf.localsfs.vss.ignoress = true;

% config for real array
conf.dimension = '2D';
conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.number = 56;
conf.secondary_sources.size = 6;
conf.secondary_sources.center = [0, 2, 0];
conf.driving_functions = 'default';
conf.xref = conf.localsfs.vss.center;

X = conf.localsfs.vss.center;
xs = [0, -1.0, 0];  % propagation direction of plane wave
src = 'pw'; 
xrange = [-conf.localsfs.vss.size/2, conf.localsfs.vss.size/2];
yrange = [-conf.localsfs.vss.size/2, conf.localsfs.vss.size/2];
zrange = 0;
f = 1000;

%%
[P_orig, x, y, z] = sound_field_mono(xrange,yrange,zrange,[xs, 1 0 0, 1],src,1,f,conf);
plot_sound_field(P_orig,x,y,z,[xs, 1 0 0, 1],conf);

conf.plot.usedb = true;
conf.usenormalisation = false;
conf.plot.colormap = 'jet';
conf.plot.caxis = [-15 1];

for r=[1/10, 1.0, 10]
% for r=[1/10, 1/3, 0.5, 2/3, 1.0, 1.5, 2.0, 3.0, 10]
  conf.localsfs.vss.logratio = r;  
  x0 = virtual_secondary_source_positions([],xs,src,conf);  
  D = driving_function_mono_wfs(x0,xs,src,f,conf);
  
  x0(:,7) = x0(:,7)./conf.localsfs.vss.number;
  [P, x, y, z] = sound_field_mono(xrange,yrange,zrange,x0,'ls',D,f,conf);
  P = P./P_orig;
  plot_sound_field(P,x,y,z,x0,conf);
  title(['ratio between shortest and longest distance: ', num2str(r)]);  
  
  x0 = secondary_source_amplitudecorrection(x0);
  P = sound_field_mono(xrange,yrange,zrange,x0,'ls',D,f,conf);
  P = P./P_orig;
  plot_sound_field(P,x,y,z,x0,conf);  
  title(['ratio between shortest and longest distance: ', num2str(r)]);  
end
