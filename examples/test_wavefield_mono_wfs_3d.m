%% Test File for 3d plot
clc
clear all
close all
%% properties of SFS_conf
conf = SFS_config;
conf.zreferenceaxis = 'y'; % set a plot reference axis
conf.useplot = 1; % plot the results = 1 ; otherwise = 0
conf.array = 'spherical'; % array type
conf.usetapwin = 0; % do not use tapering window, it's isn't needed in the 3D case
conf.plot.loudspeakers = 0; % do not plot loudspeakers in the 3D case, because it's a mess ;)
conf.usehpre = 1; % use preequalization filter
conf.number_of_points_on_sphere = 81^2; % number of points on the sphere, if spherical array is choosen
conf.debug = 1; % debug=1 allows to plot results of different evaluation steps
conf.frame = 0; % set a frame to show the wavefield in the time domain
conf.xref = [0 0 0]; % ps: 'listener position' ; pw:  place where the wavefield is scaled to one
%% properties of desired wavefield
xs = [0 1 0]; % position of point source / inicidence angle of plane wave
r = 1.5; % radius of the sphere
L = 2*r; % diameter of the sphere
src = 'fs'; % select source type pw/ps/fs
f = 1000; % evaluation frequency
conf.grid = 'HRTFgrid-';
%% plot properties
scale_axis_1 = [0 0];
scale_axis_2 = [-2 2];
scale_axis_3 = [0.5 0.5];

%% calculate frequency and time domain wavefields
wave_field_mono_wfs_3d(scale_axis_2,scale_axis_2,scale_axis_1,xs,src,f,L,conf);
% wave_field_mono_wfs_3d(scale_axis_2,scale_axis_3,scale_axis_2,xs,src,f,L,conf);
% wave_field_mono_wfs_3d(scale_axis_1,scale_axis_2,scale_axis_2,xs,src,f,L,conf);

wave_field_imp_wfs_3d(scale_axis_2,scale_axis_2,scale_axis_1,xs,src,L,conf);
grid on
% wave_field_imp_wfs_3d(scale_axis_2,scale_axis_3,scale_axis_2,xs,src,L,conf);
% grid on
% wave_field_imp_wfs_3d(scale_axis_1,scale_axis_2,scale_axis_2,xs,src,L,conf);
% grid on
%%

