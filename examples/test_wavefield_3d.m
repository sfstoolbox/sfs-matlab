%% test wave field in frequency domain
clc 
clear all
close all
%% properties
conf = SFS_config;
L = 3;
f = 1000;
conf.useplot = 1;


conf.array = 'circle';
conf.dx0 = 0.15;
X = [-2,2];
Y = [-2,2];
conf.xref = [0,0,0];
% Plane wave
src = 'pw';
xs = [0.5,1];
%% for 2.5d WFS circular array
% in frequency domain
[x,y,P_wfs25d_circular_pw,ls_activity] = wave_field_mono_wfs_25d(X,Y,xs,src,f,L,conf);
title('WFS 2.5D circular array, plane wave [0.5,1], mono-frequent');

% spatio-temporal impulse response
conf.frame = 1;
[x,y,p_wfs25d_circular_pw,ls_activity] = wave_field_imp_wfs_25d(X,Y,xs,src,L,conf);
title('WFS 2.5D circular array, plane wave [0.5,1], impulse response');
%% for 3d WFS spherical array
% conf = SFS_config;
conf.weights = [];
% diameter of array
conf.L = 3;
% desired points on the sphere
conf.number_of_points_on_sphere = 900;
% structure of secondary sources
conf.array = 'spherical';
% position of virtual source
% xs = R*[cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
xs = [0.5,1,0];

conf.usetapwin = 0;
conf.xref = [0,0,0];
conf.plot.loudspeakers = 0;

% in frequency domain
[x,y,P,win] = wave_field_mono_wfs_3d(X,Y,xs,src,f,conf.L,conf);
title('WFS 3D spherical array, plane wave [0.5,1,0], mono-frequent');
% spatio-temporal impulse response
conf.frame = [0];
[x,y,p,dds] = wave_field_imp_wfs_3d(X,Y,xs,src,L,conf);
title('WFS 3D circular array, plane wave [0.5,1,0], impulse response');
%%