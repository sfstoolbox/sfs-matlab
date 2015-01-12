%% initialize
clear variables;


%% Parameters
conf = SFS_config_example;
conf.plot.useplot = false;
conf.showprogress = true;
conf.resolution = 400;
conf.plot.loudspeakers = true;
conf.plot.realloudspeakers = false;
conf.usetapwin = true;


%% === Circular secondary sources ===
% config for real loudspeaker array
conf.dimension = '2D';
conf.secondary_sources.geometry = 'circular';
conf.secondary_sources.number = 56;
conf.secondary_sources.size = 2;
conf.secondary_sources.center = [0, 0, 0];
conf.driving_functions = 'default';
conf.xref = conf.secondary_sources.center;
% listening area
X = [0 0 0];
xs = [1.0, -1.0, 0];  % propagation direction of plane wave
src = 'pw';
xrange = [-1, 1];
yrange = [-1, 1];
zrange = 0;
f = 7000;
% --- Circular virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0, 0];
conf.localsfs.vss.geometry = 'circular';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)
% --- Linear virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0.2, 0];
conf.localsfs.vss.geometry = 'linear';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)
% --- Box shaped virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0, 0];
conf.localsfs.vss.geometry = 'box';
conf.localsfs.vss.number = 4*56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)


%% ===  Linear secondary sources ===
% config for real loudspeaker array
conf.dimension = '2D';
conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.number = 56;
conf.secondary_sources.size = 2;
conf.secondary_sources.center = [0, 1, 0];
conf.driving_functions = 'default';
% listening area
X = [0 0 0];
conf.xref = X;
xs = [1.0, -1.0, 0];  % propagation direction of plane wave
src = 'pw';
xrange = [-1, 1];
yrange = [-1, 1];
zrange = 0;
f = 7000;
% --- Circular virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0, 0];
conf.localsfs.vss.geometry = 'circular';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)
% --- Linear virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0.2, 0];
conf.localsfs.vss.geometry = 'linear';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)
% --- Box shaped virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0, 0];
conf.localsfs.vss.geometry = 'box';
conf.localsfs.vss.number = 4*56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)


%% ===  Box shaped secondary sources ===
% config for real loudspeaker array
conf.dimension = '2D';
conf.secondary_sources.geometry = 'box';
conf.secondary_sources.number = 4*56;
conf.secondary_sources.size = 2;
conf.secondary_sources.center = [0, 0, 0];
conf.driving_functions = 'default';
% listening area
X = [0 0 0];
conf.xref = X;
xs = [1.0, -1.0, 0];  % propagation direction of plane wave
src = 'pw';
xrange = [-1, 1];
yrange = [-1, 1];
zrange = 0;
f = 7000;
% --- Circular virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0, 0];
conf.localsfs.vss.geometry = 'circular';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)
% --- Linear virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0.2, 0];
conf.localsfs.vss.geometry = 'linear';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)
% --- Box shaped virtual secondary sources ---
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = false;
conf.localsfs.vss.size = 0.4;
conf.localsfs.vss.center = [0, 0, 0];
conf.localsfs.vss.geometry = 'box';
conf.localsfs.vss.number = 4*56;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
sound_field_mono_localwfs(xrange, yrange, zrange,xs,src,f,conf)
