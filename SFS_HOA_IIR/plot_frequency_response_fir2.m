
figure; % Plot...
%          ... (1) desired magnitude response vs. calculated,
%          ... (2) desired phase response vs. calculated,
%          ... (3) impulse response
%          ... (4) group delay of FIR-Filter
         
%=========================(1)=============================================
subplot(2,1,1)
plot(F, 20*log10( abs( X ) ), 'b-')
hold on
plot(config.f2,20*log10( abs( E(:,config.cLS) ) ), 'r--')
xlabel('f / Hz');
ylabel('magnitude / dB');
title('magnitude response')
grid on
%=========================(2)=============================================
subplot(2,1,2)
Q = unwrap( angle( E(:,config.cLS) ) );
plot(config.k2/(2*pi)*config.c,Q,'r');
hold on
W = unwrap((angle( X )*180./pi ));
plot(F,W,'b');
xlabel('f / Hz');
ylabel('angle / °');
title('phase response');
axis([0,7000,-350,250]);
grid on
% %=========================(3)=============================================
% subplot(2,2,3)
% plot(OUT(:,config.cLS))
% title('impulse response');
% xlabel('samples');
% ylabel('amplitude');
% %%=========================(4)==============================================
% subplot(2,2,4)
% plot(F2,20*log10(Gd))
% xlabel('f / Hz');
% ylabel('amplitude / dB');
% title('group delay of FIR Filter');