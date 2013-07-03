%TEST_WAVEFIELD_3D_WFS is a test file. No input parameters are necessary.
% Within the file different properties can be adjusted, e.g. computation of
% a plane wave ('pw'), point source ('ps') or focused source ('fs'). The 
% listener position can be adjusted by conf.xref. The used grid for 3d WFS
% is a spherical one defined by 'MaximumEnergyPoints' or the measurement
% points of the HRIRs (measured by FABIAN). 
% The frequency of the wavefield in the frequency domain can be adjusted to
% arbitrary values by editing f.

%% ===== Configuration ===================================================
clc
clear all
close all
% properties of SFS_conf
conf = SFS_config_example;
conf.plot.useplot = 1; % plot the results = 1 ; otherwise = 0
conf.plot.loudspeakers = 0; % do not plot loudspeakers in the 3D case, because it's a mess ;)
r = 1.6883; % radius of the sphere
conf.secondary_sources.geometry = 'spherical'; % array type
conf.secondary_sources.number = 81^2; 
conf.secondary_sources.grid = 'equally_spaced_points';
conf.secondary_sources.size = 2*r;
conf.usetapwin = 0; % do not use tapering window, it's not needed in the 3D case
conf.wfs.usehpre = 1; % use preequalization filter
conf.wfs.hprefhigh = 12000;
conf.debug = 1; % debug=1 allows to plot results of different evaluation steps
conf.xref = [0 0 0]; % ps/fs: 'listener position' ; pw:  place where the wavefield is scaled to one
% properties of desired wavefield
xs = [0 0.5 0]; % position of point source or focus source / inicidence angle of plane wave
src = 'fs'; % select source type pw/ps/fs
f = 1000; % evaluation frequency
conf.dimension = '3D';
% plot properties
scale_axis_1 = [0 0];
scale_axis_2 = [-2 2];
scale_axis_3 = [0.5 0.5];

%% ===== Computation ====================================================
% calculate frequency and time domain wavefields
% wave_field_mono_wfs_3d(scale_axis_2,scale_axis_2,scale_axis_1,xs,src,f,L,conf);
% title('WFS 3D spherical array r = 1.68m, focused source [0,0.5,0], mono-frequent f = 1000Hz')
% wave_field_mono_wfs_3d(scale_axis_2,scale_axis_3,scale_axis_2,xs,src,f,L,conf);
% title('WFS 3D spherical array r = 1.5m, plane wave [0,1,0], mono-frequent f = 1000Hz')
% wave_field_mono_wfs_3d(scale_axis_1,scale_axis_2,scale_axis_2,xs,src,f,L,conf);
% title('WFS 3D spherical array r = 1.5m, plane wave [0,1,0], mono-frequent f = 1000Hz')

%% calculate wavefield in time domain
wave_field_imp_wfs(scale_axis_2,scale_axis_2,scale_axis_1,xs,src,250,conf);
grid on
% title('WFS 3D spherical array r = 1.5m, plane wave [0,1,0],time domain')
% wave_field_imp_wfs_3d(scale_axis_2,scale_axis_3,scale_axis_2,xs,src,L,conf);
% grid on
% wave_field_imp_wfs_3d(scale_axis_1,scale_axis_2,scale_axis_2,xs,src,L,conf);
% grid on
%%

