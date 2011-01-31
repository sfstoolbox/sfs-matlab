% plot magnitude response for chosen LS
figure;

subplot(2,2,1)
plot( config.k2/(2*pi)*config.c, 20*log10( abs(ps_b(:, 1) ) ) );
xlabel('f / Hz');
ylabel('magnitude / dB');
title('LS bei \alpha_{0} = 0');
grid on

subplot(2,2,2)
plot( config.k2/(2*pi)*config.c, 20*log10( abs(ps_b(:, 14) ) ) );
xlabel('f / Hz');
ylabel('magnitude / dB');
title('LS bei \alpha_{0} = \pi/2');
grid on

subplot(2,2,3)
plot( config.k2/(2*pi)*config.c, 20*log10( abs(ps_b(:, 28) ) ) );
xlabel('f / Hz');
ylabel('magnitude / dB');
title('LS bei \alpha_{0} = \pi');
grid on

subplot(2,2,4)
plot( config.k2/(2*pi)*config.c, 20*log10( abs(ps_b(:, 42) ) ) );
xlabel('f / Hz');
ylabel('magnitude / dB');
title('LS bei \alpha_{0} = 3\pi/2');
grid on

%plot phase response for chosen LS
figure(2)

subplot(2,2,1)
Q = unwrap( angle( ps_b(:,1) ) );
plot(config.k2/(2*pi)*config.c, Q);
xlabel('f / Hz');
ylabel('angle / rad');
title('LS bei \alpha_{0} = 0');
grid on

subplot(2,2,2)
Q = unwrap( angle( ps_b(:,14) ) );
plot(config.k2/(2*pi)*config.c, Q);
xlabel('f / Hz');
ylabel('angle / rad');
title('LS bei \alpha_{0} = \pi/2');
grid on

subplot(2,2,3)
Q = unwrap( angle( ps_b(:,28) ) );
plot(config.k2/(2*pi)*config.c, Q);
xlabel('f / Hz');
ylabel('angle / rad');
title('LS bei \alpha_{0} = \pi');
grid on

subplot(2,2,4)
Q = unwrap( angle( ps_b(:,42) ) );
plot(config.k2/(2*pi)*config.c, Q);
xlabel('f / Hz');
ylabel('angle / rad');
title('LS bei \alpha_{0} = 3\pi/2');
grid on