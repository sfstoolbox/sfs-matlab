%% test wave field in frequency domain
clc 
clear all
close all

%% overall properties (2D/3D)

conf = SFS_config;      % load config structure
L = 2;                  % diameter of the circle/sphere
f = 1000;               % frequency for which the driving function will be calculated in frequency domain 
conf.useplot = 1;       % useplots on/off = 1/0
conf.zreferenceaxis='x';% define reference axes for all possible views
conf.dx0 = 0.15;        % distance between two secondary sources
X = [-2,2];             % begin-/endpoint of x-axis
Y = [-2,2];             % begin-/endpoint of y-axis
Z = [0,0];              % begin-/endpoint of z-axis
conf.xref = [0,1,0];    % reference point
% Plane wave
src = 'pw';             % virtual source type
xs = [0,2,0];           % coordinates of virtual source

%% optional properties

conf.debug = 1;         % debug mode on/off = 1/0;

%% for 2.5d WFS linear array
% conf.array = 'linear';
% conf.usetapwin = 0;
% conf.plot.loudspeakers = 1;
% % in frequency domain
% [x,y,P_wfs25d_circular_pw] = wave_field_mono_wfs_25d(X,Y,Z,xs,src,f,L,conf);
% title('WFS 2.5D linear array, plane wave [0,1], mono-frequent');
% % spatio-temporal impulse response
% conf.frame = 1;
% [x,y,p_wfs25d_circular_pw] = wave_field_imp_wfs_25d(X,Y,Z,xs,src,L,conf);
% title('WFS 2.5D linear array, plane wave [0,1], impulse response');

% %% for 2.5d WFS circular array
% 
% % necessary properties for circular case
% conf.array = 'circle';
% conf.usetapwin = 1;
% % in frequency domain
% [x,y,P_wfs25d_circular_pw,ls_activity] = wave_field_mono_wfs_25d(X,Y,Z,xs,src,f,L,conf);
% title('WFS 2.5D circular array, plane wave [0.5,1], mono-frequent');
% % spatio-temporal impulse response
% conf.frame = 1;
% [x,y,p_wfs25d_circular_pw,ls_activity] = wave_field_imp_wfs_25d(X,Y,Z,xs,src,L,conf);
% title('WFS 2.5D circular array, plane wave [0.5,1], impulse response');
% 
%% 3d WFS spherical array

% necessary properties for spherical case
conf.number_of_points_on_sphere = 900; % desired points on the sphere
conf.usetapwin = 0;                    % don't use tapering window, because it's not needed in the 3d spherical case
conf.plot.loudspeakers = 0;            % don't plot loudspeakers
conf.array = 'spherical';              % structure of secondary sources
conf.usehpre = 0;                      % use preequalization filter for 3D WFS
% in frequency domain
[x,y,P,win] = wave_field_mono_wfs_3d(X,Y,Z,xs,src,f,L,conf);
title('WFS 3D spherical array, plane wave [0,1,0], mono-frequent');
% spatio-temporal impulse response
conf.frame = [1];
[x,y,p,dds] = wave_field_imp_wfs_3d(X,Y,Z,xs,src,L,conf);
grid on
title('WFS 3D spherical array, plane wave [0,1,0], time domain at frame 0');