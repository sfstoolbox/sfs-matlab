%plot magnitude response for chosen LS
figure;

subplot(2,2,1)
plot( config.k2/(2*pi)*config.c, 20*log10( abs( pw_b_r(:, 1) ) ) );
xlabel('f / Hz');
ylabel('magnitude / dB');
title('LS bei \alpha_{0} = 0')
grid on

subplot(2,2,2)
plot( config.k2/(2*pi)*config.c, 20*log10( abs( pw_b_r(:, 14) ) ) );
xlabel('f / Hz');
ylabel('magnitude / dB');
title('LS bei \alpha_{0} = \pi/2')
grid on

subplot(2,2,3)
plot( config.k2/(2*pi)*config.c, 20*log10( abs( pw_b_r(:, 28) ) ) );
xlabel('f / Hz');
ylabel('magnitude / dB');
title('LS bei \alpha_{0}=\pi')
grid on

subplot(2,2,4)
plot( config.k2/(2*pi)*config.c, 20*log10( abs( pw_b_r(:, 42) ) ) );
xlabel('f / Hz');
ylabel('magnitude / dB');
title('LS bei \alpha_{0}=3\pi/2')
grid on

%plot phase response for chosen LS
figure;
axis([0,7000,0,360])
subplot(2,2,1)
Q = unwrap( angle( pw_b_r(:,1) ) );
plot(config.k2/(2*pi)*config.c, Q);
xlabel('f / Hz');
ylabel('angle / rad');
title('LS bei \alpha_{0} = 0');
axis([0,7000,-100,250])
grid on

subplot(2,2,2)
Q = unwrap( angle( pw_b_r(:,14) ) );
plot(config.k2/(2*pi)*config.c, Q);
xlabel('f / Hz');
ylabel('angle / rad');
title('LS bei \alpha_{0} = \pi/2');
axis([0,7000,-100,250])
grid on

subplot(2,2,3)
Q = unwrap( angle( pw_b_r(:,28) ) );
plot(config.k2/(2*pi)*config.c, Q);
xlabel('f / Hz');
ylabel('angle / rad');
title('LS bei \alpha_{0} = \pi');
axis([0,7000,-100,250])
grid on

subplot(2,2,4)
Q = unwrap( angle( pw_b_r(:,42) ) );
plot(config.k2/(2*pi)*config.c, Q);
xlabel('f / Hz');
ylabel('angle / rad');
title('LS bei \alpha_{0} = 3\pi/2');
axis([0,7000,-100,250])
grid on
