%%DEBUG_FAR_FIELD_HRIRs compares the ILDs and HRIRs of the given 3D HRIR 
% dataset of FABIAN with the range extrapolated far-field HRIRs and their
% ILDs. For ease of comparision you have to select your desired elevation
% angle with angle=''; (in degree). Comparision is possible in steps of 2°
% from -64°...88°. Therefore the FABIAN HRIR and the extrapolated HRIR dataset 
% have to be stored in arc of circles of the sphere with the file name 
% irs_theta_"angle" / irs_pw_theta_"angle" (e.g.: irs_theta_0 / irs_pw_theta_0)
% 
% see: find_arc_of_circles_for_different_elevation_angles

clear all;close all;clc

angle = '74';

load(['irs_theta_' angle '.mat'])
ild1 = interaural_level_difference(left,right);
irs_left = left;
irs_right = right;

figure 
subplot(2,2,1)
imagesc([-180:180],[0:425]/441,irs_left)
colormap(gray);
xlabel('$\varphi$ in degree','interpreter','latex')
ylabel('time (ms)')
title(['HRIRs (elevation = ' angle '°)'])

load(['irs_pw_theta_' angle '.mat']);
ild2 = interaural_level_difference(left,right);
irs_pw_left = left;
irs_pw_right = right;

subplot(2,2,2)
imagesc([-180:180],[0:425]/441,irs_pw_left)
colormap(gray);
xlabel('$\varphi$ in degree','interpreter','latex')
ylabel('time (ms)')
title(['synthetic far-field HRIRs (elevation = ' angle '°)'])

subplot(2,2,[3 4])
plot(sort(degree(apparent_azimuth),'descend'),ild1,'r')
hold on
plot(sort(degree(apparent_azimuth),'descend'),ild2,'b')
legend(['without extrapolation, elevation = ' angle '°'],['with extrapolation, elevation = ' angle '°'])
title('ILD')

xlabel('$\varphi$ in degree','interpreter','latex')
ylabel('amplitude')
grid on

% %% frequency response
% [amplitude1,phase1,f] = easyfft(irs.left(:,1),conf);
% [amplitude2,phase2,f] = easyfft(irs_pw.left(:,1),conf);
% 
% figure
% plot(f,db(amplitude1),'r','MarkerSize',10)
% hold on
% plot(f,db(amplitude2),'b','MarkerSize',10)
% grid on
