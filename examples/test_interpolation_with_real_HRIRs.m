%==========================================================================
% Test HRIR Interpolation %%get_ir(irs,phi,delta,r,X0)
clc
%==========================================================================
%% 2D Dataset (circle) --> 1° phi resolution

% irs = read_irs('QU_KEMAR_anechoic_AKGK271_3m.mat');
% irs = read_irs('QU_KEMAR_anechoic_3m.mat');
% % 211 --> 0.523598775598299
% % 212 --> 0.541052068118242
% ir = get_ir(irs,0.523598775598299,0,3,[0 0 0]);
% irs.apparent_azimuth(1,211) = 0;
% ir2 = get_ir(irs,0.523598775598299,0,3,[0 0 0]);

%% 3D Dataset (sphere) --> 7.2° phi resoultion, 7.2° theta resolution (1250 values)

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

%% 3D Datatset (sphere) --> 2° phi and theta resolution
irs = read_irs('FABIAN3d1,7.mat');
% % 405 --> phi = 216
% % 405 --> theta = 86
ir = get_ir(irs,-2.5133,1.5010,1.6520,[0 0 0]);
irs.apparent_azimuth(1,22) = 0;
irs.apparent_elevation(1,22) = 0;
irs.distance(1,22) = 0;
ir2 = get_ir(irs,-2.5133,1.5010,1.6520,[0 0 0]);

%% plots

figure
subplot(2,2,1)
plot(ir(:,1),'b')
hold on
plot(ir2(:,1),'r-.')
grid on
xlabel('samples')
ylabel('amplitude')
title('Left ear HRIRs without and with interpolation, dataset: dataset: FABIAN,3D,resolution theta 2°')
legend('without interpolation','with interpolation')

subplot(2,2,2)
plot(ir(:,2),'g')
hold on
plot(ir2(:,2),'k-.')
grid on
xlabel('samples')
ylabel('amplitude')
title('Right ear HRIRs without and with interpolation, dataset: dataset: FABIAN,3D,resolution theta 2°')
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
[amplitude_ir1r,phase_ir1r,f1r] = easyfft(ir(:,2));
[amplitude_ir2l,phase_ir2l,f2l] = easyfft(ir2(:,1));
[amplitude_ir2r,phase_ir2r,f2r] = easyfft(ir2(:,2));

figure
subplot(2,2,1)
plot(f1l,db(amplitude_ir1l),'b')
hold on
plot(f1l,db(amplitude_ir2l),'r-.')
grid on
xlabel('f/Hz')
ylabel('amplitude / db')
title('Left ear HRTFs without and with interpolation, dataset: FABIAN,3D,resolution theta 2°')
legend('without interpolation','with interpolation')

subplot(2,2,2)
plot(f2l,db(amplitude_ir1r),'g')
hold on
plot(f2l,db(amplitude_ir2r),'k-.')
grid on
xlabel('f / Hz')
ylabel('amplitude / db')
title('Right ear HRTFs without and with interpolation, dataset: dataset: FABIAN,3D,resolution theta 2°')
legend('without interpolation','with interpolation')

subplot(2,2,3)
plot(f1l,db(amplitude_ir2l)-db(amplitude_ir1l),'r')
grid on
xlabel('f / Hz')
ylabel('amplitude / db')
title('Error between measured and interpolated HRTFs (left ear)')

subplot(2,2,4)
plot(f1l,db(amplitude_ir2r)-db(amplitude_ir1r),'r')
grid on
xlabel('f / Hz')
ylabel('amplitude / db')
title('Error between measured and interpolated HRTFs (right ear)')