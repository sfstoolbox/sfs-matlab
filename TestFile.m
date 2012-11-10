%Test File
clc 
clear all
close all
%% properties which need to be set
conf = SFS_config_example;
% diameter of array
L = 1;
% structure of secondary sources
conf.array = 'spherical';
% resolution of elevation angle for spherical arrays
conf.resolution_theta = 10;

% position of virtual source
R = 1;
phi = pi/2;
theta = pi/2;
xs = R*[cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];

%% calculate position and direction of secondary sources
x0 = secondary_source_positions(L,conf);

%Activity of secondary sources
ls_activity = secondary_source_selection(x0,xs,'pw');
idx = ((ls_activity>0));

%% plot all secondary sources and active secondary sources to prove the evaluation
figure
subplot(2,2,1)
plot3(x0(:,1),x0(:,2),x0(:,3),'bx')
hold on
plot3(x0(idx,1),x0(idx,2),x0(idx,3),'rx')
hold on
plot3(xs(1),xs(2),xs(3),'gx','Linewidth',10)
 xlabel('x->[m]')
 ylabel('y->[m]')
 zlabel('z->[m]')
 grid on
 axis equal
 legend('non-active secondary sources','active secondary sources','incidence angle of pw')
 title('Secondary Sources of 3D WFS')
 
subplot(2,2,4)
plot3(x0(idx,1),x0(idx,2),x0(idx,3),'rx')
 hold on
plot3(xs(1),xs(2),xs(3),'gx','Linewidth',10)
 xlabel('x->[m]')
 ylabel('y->[m]')
 zlabel('z->[m]')
 grid on
 axis equal
 title('Active Secondary Sources of 3D WFS for desired \theta,\phi')
 
subplot(2,2,3)
plot3(x0(:,1),x0(:,2),x0(:,3),'bx')
 xlabel('x->[m]')
 ylabel('y->[m]')
 zlabel('z->[m]')
 grid on
 axis equal
 title('All Secondary Sources of 3D WFS')
 %%