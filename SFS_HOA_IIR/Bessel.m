function [w,config] = Bessel(config)
%==========================================================================
% show hankelfunction
x = 0:0.001:20;

o = sphbesselh(0, 2, x);
p = sphbesselh(1, 2, x);
q = sphbesselh(2, 2, x);
z = hankel_rekursiv_show(config);

figure;
subplot(2,2,1)
plot(x,real(o));
hold on
plot(x,real(p),'r');
hold on
plot(x,real(q),'g');
%title('real part')
grid on;
axis([0 20 -0.5 1])


subplot(2,2,2)
plot(x,imag(o));
hold on
plot(x,imag(p),'r');
hold on
plot(x,imag(q),'g');
grid on;
axis([0 20 -1 3])
%title('imaginary part')

subplot(2,2,3)
plot(x,20*log10(abs(o)))
hold on
plot(x,20*log10(abs(p)),'r')
hold on
plot(x,20*log10(abs(q)),'g')
%title('absolute value')
grid on;
axis([0 20 -30 30])


subplot(2,2,4)
plot(x,20*log10(abs(1./o)))
hold on
plot(x,20*log10(abs(1./p)),'r')
hold on
plot(x,20*log10(abs(1./q)),'g')
%title('absolute value')
grid on;
axis([0 20 -30 30])

% figure;
% subplot(2,2,1)
% plot(x,real(q))
% hold on
% plot(real(z(config.p+1,:)),'r:','LineWidth',2)
% %==========================================================================
% % non recursive calculation
% y = sphbesselh(config.p, 2, config.k2.*config.R );

% % recursive calculation of Hankelfunction
% z = hankel_rekursiv(config);
% z = hankel_rekursiv_real_imag(config);
% 
% figure;
% subplot(2,2,1)
% plot(config.f2,20*log10(abs(y)))
% hold on
% plot(config.f2,20*log10(abs(z(config.p+1,:))),'r:','LineWidth',2);
% title('absolute value of hankelfunction 2nd kind non rec (blue) and rek (red)')
% axis([0 20000 -30 30])
% 
% subplot(2,2,2)
% plot(config.f2,real(y));
% hold on
% plot(config.f2,real(z(config.p+1,:)),'r:','LineWidth',2);
% title('real part of hankelfunction 2nd kind non rec (blue) and rek (red)')
% 
% subplot(2,2,3)
% plot(config.f2,imag(y));
% hold on
% plot(config.f2,imag(z(config.p+1,:)),'r:','LineWidth',2);
% title('imaginary part of hankelfunction 2nd kind non rec (blue) and rek (red)')
% 

w = 1;
end



