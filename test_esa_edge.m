SFS_start;

conf = SFS_config;
conf.secondary_sources.size = 3;
conf.secondary_sources.number = 512;
conf.secondary_sources.geometry = 'edge';
conf.secondary_sources.center = [0,0,0];
conf.secondary_sources.alpha = [-1/2*pi, pi];

conf.plot.realloudspeakers = false;
conf.plot.useplot = true;
conf.plot.usenormalisation = true;
conf.plot.normalisation = 'center';

f = 1000;
xs = [3,-1,0];

x0 = secondary_source_positions(conf);

D = driving_function_mono_esa_edge(x0,xs,'ls',f,conf);

X = [-3,0];
Y = [-3,0];
Z = 0;

[P,x,y,z] = sound_field_mono(X,Y,Z,x0,'ls',D,f,conf);
P_gt = sound_field_mono_line_source(X,Y,Z,xs,f,conf);

% figure, plot(x, db([P, P_gt]))