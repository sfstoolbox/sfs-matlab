%% Test File for 3d plot
clc
clear all
close all
%%
conf = SFS_config;
conf.zreferenceaxis = 'y';
conf.useplot = 1;
conf.array = 'spherical';
conf.usetapwin = 0;
conf.plot.loudspeakers = 0;
conf.usehpre = 1;
conf.number_of_points_on_sphere = 20^2;
conf.debug = 1;
[x,y,P,win] = wave_field_mono_wfs_3d([-3 3],[-0.15 3],[0 0],[0 1],'pw',1000,3,conf);
% [x,y,P,win] = wave_field_mono_wfs_3d([-3 3],[0.5 0.5],[-3 3],[0 1 0],'pw',1000,3,conf);
% [x,y,P,win] = wave_field_mono_wfs_3d([0 0],[-0.15 3],[-3 3],[0 1 0],'pw',1000,3,conf);
%%

