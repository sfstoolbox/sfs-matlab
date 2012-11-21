%% Test File to prove secondary source position for spherical array and
%  active LS for 3D WFS
clc 
clear all
close all

%% properties 
conf = SFS_config;
% diameter of array
conf.L = 1;
% desired points on the sphere
conf.number_of_points_on_sphere = 900;
% structure of secondary sources
conf.array = 'spherical';
% position of virtual source
R = 1;
phi = 0;
theta = pi/4;
xs = R*[cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];

%% calculate secondary source positions and directions
x0 = secondary_source_positions(conf.L,conf);

%% Calculate active secondary sources
x0 = x0(:,1:6);
x0 = secondary_source_selection(x0,xs,'pw',conf.xref);

%% plot all secondary sources and active secondary sources to prove the evaluation
figure
plot3(x0(:,1),x0(:,2),x0(:,3),'bx')
% hold on
% plot3(x0(idx,1),x0(idx,2),x0(idx,3),'rx')
hold on
plot3(xs(1),xs(2),xs(3),'gx','Linewidth',10)
 xlabel('x->[m]')
 ylabel('y->[m]')
 zlabel('z->[m]')
 grid on
 axis([-conf.L conf.L -conf.L conf.L -conf.L conf.L])
 axis equal
 legend('active secondary sources','incidence angle of pw')
 title('Secondary Sources of 3D WFS')
 
% subplot(2,2,4)
% plot3(x0(idx,1),x0(idx,2),x0(idx,3),'rx')
%  hold on
% plot3(xs(1),xs(2),xs(3),'gx','Linewidth',10)
%  xlabel('x->[m]')
%  ylabel('y->[m]')
%  zlabel('z->[m]')
%  grid on
%  axis equal
%  title('Active Secondary Sources of 3D WFS for desired \theta,\phi')
%  
% subplot(2,2,3)
% plot3(x0(:,1),x0(:,2),x0(:,3),'bx')
%  xlabel('x->[m]')
%  ylabel('y->[m]')
%  zlabel('z->[m]')
%  grid on
%  axis equal
%  title('All Secondary Sources of 3D WFS')
