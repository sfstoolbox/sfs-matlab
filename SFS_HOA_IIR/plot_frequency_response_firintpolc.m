%==========================================================================
figure;   % Plot...
            % ... (1) desired magnitude response vs. calculated,
            % ... (2) desired phase response vs. calculated,
            % ... (3) impulse response
            % ... (4) group delay of FIR-Filter
%=========================(1)==============================================           
subplot(2,1,1)
plot(F, 20*log10( abs( X ) ),'b--')
hold on
plot(config.f2, 20*log10( abs( E(:,config.cLS) ) ),'r-.','LineWidth',2 )
xlabel('f / Hz');
ylabel('magnitude / dB');
title('magnitude response')
grid on
%=========================(2)==============================================
subplot(2,1,2)
Q = unwrap( angle( E(:,config.cLS) ) );
 plot(config.f2,Q ,'r-.','LineWidth',2);
 hold on
 W = unwrap((angle( (X) )));
 plot(F,W,'b--');
 xlabel('f / Hz');
 ylabel('angle / rad');
 title('phase response');
 grid on
 %axis([0,20000,-1000,600])
% %=========================(3)==============================================
% subplot(2,2,3)
% plot( OUT(:,config.cLS) )
% title('impulse response');
% %hold on
% %plot(T1,H1,'r');
% %=========================(4)==============================================
% subplot(2,2,4)
% plot(F2,Gd,'r')
% xlabel('f / Hz');
% ylabel('amplitude / dB');
% title('group delay of FIR Filter');
% %==========================================================================