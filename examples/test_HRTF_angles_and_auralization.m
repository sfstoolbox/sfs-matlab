%% check angles of HRTFs and auralize HRTFs
clear all;close all; clc;

%% CHECK ANGLES OF THE HRTF DATASET
%% load HRTF dataset
% load('az0.mat');
% phi_d = source_azimuth(:,:)';
% theta_d = source_elevation';
% r_original = distance';

%% angles in rad
% phi_rad = rad(phi_d);
% theta_rad = rad(theta_d);

%% angles in rad and corrected
% phi_corrected = correct_azimuth(rad(phi_d));
% theta_corrected = correct_elevation(rad(theta_d));

%% compare the different angles
% % compare = [phi_d ;phi_rad;phi_corrected ;zeros(1,length(distance));...
%             theta_d;theta_rad;theta_corrected;];
%         
%% in cartesian coordinates
% [x_original,y_original,z_original] = sph2cart(phi_d,theta_d,r_original);
% [x_rad,y_rad,z_rad] = sph2cart(phi_rad,theta_rad,r_original);
% [x_corrected,y_corrected,z_corrected] = sph2cart(phi_corrected,theta_corrected,r_original);

% compare_places = [x_original;y_original;z_original;zeros(1,length(r_original)); ...
%                   x_rad;y_rad;z_rad;zeros(1,length(r_original));...
%                   x_corrected;y_corrected;z_corrected];

%% AURALIZATION
%% create sinus for auralization with HRTFs
t = 0:1/44100:10;
sinus = 2*sin(2*pi*440*t);
wavwrite([sinus' sinus'],44100,32,'Sinus440Hz');

%% load dataset and apply properties
conf = SFS_config;
conf.array = 'spherical';
conf.usehpre = 1;
conf.fs = 44100;
load('FABIAN_3D_anechoic_~1.6.mat');

%% choose the desired point
point_phi = correct_azimuth(rad(0));
point_theta = correct_elevation(rad(0));

%% get radius because r is various
idx = findrows(...
    [irs.apparent_azimuth' irs.apparent_elevation'],...
    [point_phi,point_theta]);
R = irs.distance(1,idx);

%% auralization
ir = get_ir(irs,point_phi,point_theta,R,[0 0 0]);
x = wavread('Sinus440Hz');
outsig = auralize_ir(ir,x,1,conf);
wavwrite(outsig,44100,32,'HRTF_270');

%% plot ILD
ild1 = interaural_level_difference(ir(:,1),ir(:,2));
figure
plot(ild1,'rx','MarkerSize',10)
hold on
grid on