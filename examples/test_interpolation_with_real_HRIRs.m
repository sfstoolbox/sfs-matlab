%==========================================================================
% Test HRIR Interpolation %%get_ir(irs,phi,delta,r,X0)
clc
clear all
close all
%==========================================================================
%% 2D Dataset (circle) --> 1° phi resolution

% irs = read_irs('QU_KEMAR_anechoic_AKGK271_3m.mat');
% irs = read_irs('QU_KEMAR_anechoic_3m.mat');
% % 211 --> 0.523598775598299
% % 212 --> 0.541052068118242
% ir = get_ir(irs,0.523598775598299,0,3,[0 0 0]);
% irs.apparent_azimuth(1,211) = 0;
% ir2 = get_ir(irs,0.523598775598299,0,3,[0 0 0]);

%% 3D Dataset (sphere) --> 7.2° resoultion, 7.2° theta resolution (1250 values)

% irs = read_irs('CIPIC_KEMAR_large.mat');
% % 874 --> phi = 0.7854
% % 874 --> theta = 1.4726
% ir = get_ir(irs,0.7854,1.4726,1,[0 0 0]);
% irs.apparent_azimuth(1,874) = 0;
% irs.apparent_elevation(1,874) = 0;
% ir2 = get_ir(irs,0.7854,1.4726,1,[0 0 0]);

%% 3D Dataset (sphere) --> 10° phi resolution , 10° theta resolution (705 values)

% irs = read_irs('MIT_KEMAR_anechoic_1.7m.mat');
% % 405 --> phi = 0.4189
% % 405 --> theta = 0.5236 
% ir = get_ir(irs,0.4189,0.5236,1.7,[0 0 0]);
% irs.apparent_azimuth(1,405) = 0;
% irs.apparent_elevation(1,405) = 0;
% ir2 = get_ir(irs,0.4189,0.5236,1.7,[0 0 0]);

%% 3D Datatset (sphere) --> r~1.6m

% load HRIR data set
basepath = which('test_interpolation_with_real_HRIRs');
basepath = basepath(1:end-36);
addpath([basepath,'3D_HRIR_DataSet_1,6']);
irs = read_irs('FABIAN_3D_anechoic_~1.6.mat');
% get random HRIR of data set
random_ir = randi(length(irs.left),1,1);
% random_ir = 5886;% --> (0,0)
% calculate the cartesian coordinates of the place of the random chosen HRIR 
[x1,y1,z1] = sph2cart(irs.apparent_azimuth(1,random_ir),irs.apparent_elevation(1,random_ir),irs.distance(1,random_ir));

% search for the HRIR in the data set
ir = get_ir(irs,irs.apparent_azimuth(1,random_ir),irs.apparent_elevation(1,random_ir),irs.distance(1,random_ir),[0 0 0]);

% defining help variables to store the measurement point of the random chosen HRIR
help1 = irs.apparent_azimuth(1,random_ir);
help2 = irs.apparent_elevation(1,random_ir);
help3 = irs.distance(1,random_ir);

% set the measurement point to an arbitrary negative value
irs.apparent_azimuth(1,random_ir) = -1;
irs.apparent_elevation(1,random_ir) = -1;
irs.distance(1,random_ir) = -1;

% try to get the HRIR at the random chosen point again --> now
% interpolation is necessary
ir2 = get_ir(irs,help1,help2,help3,[0 0 0]);

% if debug = 1, plots will be shown
debug = 1;

%% plots
if debug
    
figure
subplot(2,2,1)
plot(ir(:,1),'b')
hold on
plot(ir2(:,1),'r-.')
grid on
xlabel('samples')
ylabel('amplitude')
title('Left ear HRIRs without and with interpolation, dataset: dataset: FABIAN,3D r~1.6m')
legend('without interpolation','with interpolation')

subplot(2,2,2)
plot(ir(:,2),'g')
hold on
plot(ir2(:,2),'k-.')
grid on
xlabel('samples')
ylabel('amplitude')
title('Right ear HRIRs without and with interpolation, dataset: dataset: FABIAN,3D,r~1.6m')
legend('without interpolation','with interpolation')

subplot(2,2,3)
plot(abs(ir2(:,1)-ir(:,1)),'r')
grid on
xlabel('samples')
ylabel('amplitude')
title('Error between measured and interpolated HRIR (left ear)')

subplot(2,2,4)
plot(abs(ir2(:,2)-ir(:,2)),'r')
grid on
xlabel('samples')
ylabel('amplitude')
title('Error between measured and interpolated HRIR (right ear)')

%% HRTF
[amplitude_ir1l,phase_ir1l,f1l] = easyfft(ir(:,1));
[amplitude_ir1r,phase_ir1r,~] = easyfft(ir(:,2));
[amplitude_ir2l,phase_ir2l,~] = easyfft(ir2(:,1));
[amplitude_ir2r,phase_ir2r,~] = easyfft(ir2(:,2));

figure
subplot(3,2,1)
plot(f1l,db(amplitude_ir1l),'b')
hold on
plot(f1l,db(amplitude_ir2l),'r-.')
grid on
xlabel('f/Hz')
ylabel('amplitude / db')
title('Left ear HRTFs without and with interpolation, dataset: FABIAN,3D,r~1.6m')
legend('without interpolation','with interpolation')

% HRTF from original data set
% hold on
% load('az0.mat')
% [amplitude_original1,phase_original1,f1l] = easyfft(left(:,random_ir));
% plot(f1l,db(amplitude_original1),'m-.')
% [amplitude_original2,phase_original2,~] = easyfft(right(:,random_ir));


subplot(3,2,2)
plot(f1l,db(amplitude_ir1r),'g')
hold on
plot(f1l,db(amplitude_ir2r),'k-.')
grid on
xlabel('f / Hz')
ylabel('amplitude / db')
title('Right ear HRTFs without and with interpolation, dataset: dataset: FABIAN,3D,r~1.6m')
legend('without interpolation','with interpolation')

% hold on
% plot(f1l,db(amplitude_original2),'m-.')

subplot(3,2,3)
plot(f1l,db(amplitude_ir2l)-db(amplitude_ir1l),'r')
grid on
xlabel('f / Hz')
ylabel('amplitude / db')
title('Error between measured and interpolated HRTFs (left ear)')

subplot(3,2,4)
plot(f1l,db(amplitude_ir2r)-db(amplitude_ir1r),'r')
grid on
xlabel('f / Hz')
ylabel('amplitude / db')
title('Error between measured and interpolated HRTFs (right ear)')
%% Proof the nearest points
load('NearestPoints.mat')
subplot(3,2,5)
[x2,y2,z2] = sph2cart(irs.apparent_azimuth(1,idx(1)),irs.apparent_elevation(1,idx(1)),irs.distance(1,idx(1)));
[x3,y3,z3] = sph2cart(irs.apparent_azimuth(1,idx(2)),irs.apparent_elevation(1,idx(2)),irs.distance(1,idx(2)));
[x4,y4,z4] = sph2cart(irs.apparent_azimuth(1,idx(3)),irs.apparent_elevation(1,idx(3)),irs.distance(1,idx(3)));
plot3(x1,y1,z1,'rx','MarkerSize',10)
hold on
plot3(x2,y2,z2,'bx','MarkerSize',10)
hold on
plot3(x3,y3,z3,'bx','MarkerSize',10)
hold on
plot3(x4,y4,z4,'bx','MarkerSize',10)
hold on
grid on 
N =20;
phi = linspace(0,(1-1/N)*2*pi,N);
theta = linspace(-pi,(1-1/N)*pi,N);
[PHI,THETA] = meshgrid(phi,theta);
r=1.6;
x = r.*cos(PHI).*sin(THETA);
y = r.*sin(PHI).*sin(THETA);
z = r.*cos(THETA);

a=mesh(x,y,z,zeros(size(z)));
colormap(gray)
daspect([1 1 1])
title('nearest neighbours (blue) of desired interpolation point (red)')
xlabel('x --> [m]')
ylabel('y --> [m]')
zlabel('z --> [m]')

%% considering ILD
subplot(3,2,6)
plot(f1l,db(amplitude_ir1l)-db(amplitude_ir1r),'r');
hold on
plot(f1l,db(amplitude_ir2l)-db(amplitude_ir2r),'b');
xlabel('f / Hz')
ylabel('amplitude / db')
title('ILD')
legend('measured HRTF','interpolated HRTF')
grid on

end