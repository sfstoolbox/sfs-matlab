%==========================================================================
% Test HRIR Interpolation %%get_ir(irs,phi,delta,r,X0)
clc
%==========================================================================
irs = read_irs('QU_KEMAR_anechoic_3m.mat');
% 211 --> 0.523598775598299
% 212 --> 0.541052068118242
% rad(30.3) = 0.5288
ir = get_ir(irs,0.523598775598299,0,1,[0 0 0]);
irs.apparent_azimuth(1,211) = 0;
ir2 = get_ir(irs,0.523598775598299,0,1,[0 0 0]);

figure
plot(ir(:,1),'b')
hold on
plot(ir2(:,1),'r')
grid on
xlabel('samples')
ylabel('amplitude')
title('Left ear HRIRs without and with interpolation')
legend('without interpolation','with interpolation')


figure
plot(ir(:,2),'g')
hold on
plot(ir2(:,2),'k')
grid on
xlabel('samples')
ylabel('amplitude')
title('Right ear HRIRs without and with interpolation')
legend('without interpolation','with interpolation')
