%% test interpolation
clc
clear all
close all

%% properties 
conf = SFS_config;
% diameter of array
L = 2;
% structure of secondary sources
conf.array = 'spherical';
conf.number_of_points_on_sphere = 81;
%  calculate position and direction of secondary sources
x0 = secondary_source_positions(L,conf);
% x0 = x0(:,1:6);
%% plot a sphere only for validation 

N =20;
phi = linspace(0,(1-1/N)*2*pi,N);
theta = linspace(-pi,(1-1/N)*pi,N);
[PHI,THETA] = meshgrid(phi,theta);
r=1;
x = r.*cos(PHI).*sin(THETA);
y = r.*sin(PHI).*sin(THETA);
z = r.*cos(THETA);

figure
a=mesh(x,y,z,zeros(size(z)));
colormap(gray)
daspect([1 1 1])
xlabel('x')
ylabel('y')
zlabel('z')
% set(a,'facecolor','none')
title('Illustration of Interpolation of HRIRs using VBAP')
hold on
%% interpolation (of HRIRs) 
% given loudpeaker positions
phi1 = pi/4;
phi2 = pi/4;
phi3 = pi/4;

theta1 = pi/4;
theta2 = pi/4;
theta3 = pi/4;

% desired angle of new virtual loudspeaker
% note: phi_min < alpha < phi_max AND theta_min < beta < theta_max
 alpha = pi/4
 beta =  pi/4

% P = 30;
% [alpha,beta,r] = cart2sph(x0(P,1),x0(P,2),x0(P,3))
% points
[ir_x1,ir_y1,ir_z1] = sph2cart(phi1,theta1,1);
[ir_x2,ir_y2,ir_z2] = sph2cart(phi2,theta2,2);
[ir_x3,ir_y3,ir_z3] = sph2cart(phi3,theta3,3);
[ir_x,ir_y,ir_z] = sph2cart(alpha,beta,2.5);

ir1 = [ir_x1,ir_y1,ir_z1];
ir2 = [ir_x2,ir_y2,ir_z2];
ir3 = [ir_x3,ir_y3,ir_z3];
ir =  [ir_x,ir_y,ir_z];

L = [ir1;ir2;ir3];
p = ir;

L
p

if rank(L) == 3
    
    g = L.'\p';
    ir = g(1)*ir1 + g(2).*ir2 + g(3).*ir3; 
    
elseif rank(L) == 2

    L = [ir1(1:2);ir2(1:2)];
    p = p(1:2);
    g = L.'\p';
    ir = g(1)*ir1 + g(2).*ir2;
    
elseif rank(L) == 1
    
    L = ir1(1);
    p = p(1);
    g = L.'\p';
    ir = g(1)*ir1;
    
end
% calculate desired ir

[phi_desired,theta_desired,r_desired] = cart2sph(ir(1,1),ir(1,2),ir(1,3))
%% plot
plot3(x0(:,1),x0(:,2),x0(:,3),'bx','MarkerSize',8)
hold on
plot3(ir_x1,ir_y1,ir_z1,'rx','LineWidth',3)
hold on
plot3(ir_x2,ir_y2,ir_z2,'bx','LineWidth',3)
hold on
plot3(ir_x3,ir_y3,ir_z3,'kx','LineWidth',3)
grid on
plot3(ir(1,1),ir(1,2),ir(1,3),'go','LineWidth',3,'MarkerSize',8)
legend('unit sphere','desired virtual secondary source grid','measured HRIR1','measured HRIR2','measured HRIR3','Interpolated HRIR at desired secondary source')