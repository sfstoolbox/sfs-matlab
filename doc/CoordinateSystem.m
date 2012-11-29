%% nomenclature --> picture
clc 
clear all 
close all
%% coordinate system
N = 25;
phi1 = pi/4;
R=linspace(0,1.2,N);
x = linspace(0,1.2,N);
y = x;
z = x;
help = zeros(1,N);
%% circle for azimuth angle
R = 1;
phi = linspace(-pi/6,(pi/2+pi/6),N);
kreis = R*[cos(phi);sin(phi)];
%% point on sphere
theta = pi/4;
phi = pi/6;
R = repmat(linspace(0,1,N),3,1);
x0 = R.*[repmat(cos(phi)*sin(theta),1,N);repmat(sin(phi)*sin(theta),1,N);repmat(cos(theta),1,N)];
%% help vectors
M=8;
R= linspace(0,1,M);
xhelp1_1 = R.*repmat(x0(1,end),1,M);
xhelp1_2 = R.*repmat(x0(2,end),1,M);
xhelp1_3 = zeros(1,M);

xhelp2_1 = repmat(xhelp1_1(end),1,M);
xhelp2_2 = repmat(xhelp1_2(end),1,M);
xhelp2_3 = theta.*linspace(0,1,M);

point = [x0(1,end) x0(2,end) x0(3,end)];
%% circle for elevation angle
R=1;
theta = linspace(pi/2+pi/16,-pi/16,N);
kreis2 = R*[cos(theta);zeros(1,N);sin(theta)];
kreis3 = rotation_matrix(pi/4)*kreis2;
kreis2 = rotation_matrix(pi/6)*kreis2;
%% angles
% theta
phi = linspace(0,pi/4,M);
theta_angle = 0.4.*[cos(phi);zeros(1,length(phi));sin(phi)];
phi = pi/6;
R = rotation_matrix(phi);
theta_angle = R*theta_angle;
% phi
phi = linspace(0,pi/6,M); 
phi_angle = 0.6*[cos(phi);sin(phi);zeros(1,length(phi))];

%% plot all stuff
legend1=[];
legend1(1) = plot3(x,help,help,'k','LineWidth',2);
hold on
legend1(2) = plot3(help,y(1,:),help,'k','LineWidth',2);
hold on
legend1(3) = plot3(help,help,z,'k','LineWidth',2);
hold on
legend1(4) = plot3(kreis(1,:),kreis(2,:),help,'k','LineWidth',2);
hold on
% plot3(kreis2(1,:),kreis2(2,:),kreis2(3,:),'k.','LineWidth',2)
% hold on
legend1(5) = plot3(kreis3(1,:),kreis3(2,:),kreis3(3,:),'k','LineWidth',2);
hold on
legend1(6) = plot3(x0(1,:),x0(2,:),x0(3,:),'k');
hold on
legend1(7) = plot3(xhelp1_1,xhelp1_2,xhelp1_3,'k.');
hold on
legend1(8) = plot3(xhelp2_1(1:end-1),xhelp2_2(1:end-1),xhelp2_3(1:end-1),'k.');
hold on
plot3(point(1),point(2),point(3),'ko','MarkerFaceColor',[0 0 0]);
text(point(1)+0.05,point(2)+0.05,point(3),'x_{0}',...
     'HorizontalAlignment','left');
hold on
legend1(9) = plot3(phi_angle(1,:),phi_angle(2,:),phi_angle(3,:),'k.-');
text(phi_angle(1,round(M/2))-0.1,phi_angle(2,round(M/2)),phi_angle(3,round(M/2)),'\phi',...
     'HorizontalAlignment','left');
hold on
legend1(10) = plot3(theta_angle(1,:),theta_angle(2,:),theta_angle(3,:),'kx-');
text(theta_angle(1,round(M/2))-0.05,theta_angle(2,round(M/2)),theta_angle(3,round(M/2))-0.1,'\theta',...
     'HorizontalAlignment','left');
% arrows and axis names 
% text(-0.05,0,1.2,'\uparrow',...
%      'HorizontalAlignment','left')
text(-0.08,0,1.2,'z',...
     'HorizontalAlignment','left')
% text(1.15,-0.05,'\rightarrow',...
%      'HorizontalAlignment','left')
text(1.17,-0.12,'x',...
     'HorizontalAlignment','center')  
text(+0.012,y(end)+0.06,0,'y',...
     'HorizontalAlignment','center')  
axis off

view([1.75,-7.4,6.5])
%% Notes
title('Used Coordinate System','FontSize',16)

l=legend(legend1([9 10]),'$-\pi$ $<$ $\Phi$ $<$ $\pi$','$\frac{\pi}{2}$ $<$ $\Theta$ $<$ $\frac{\pi}{2}$','Location','North')
set(l,'Interpreter','Latex','FontSize',16);
set(l, 'Box', 'off')
set(l, 'Color', 'none')
