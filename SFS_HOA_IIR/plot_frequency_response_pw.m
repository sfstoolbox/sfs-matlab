%plot magnitude response for chosen LS
figure;

subplot(2,2,1)
plot( config.k2/(2*pi)*config.c, 20*log10( abs( pw_b(:, 1) ) ) );
xlabel('f / Hz');
ylabel('magnitude / dB');
title('\alpha_{0} = 0')

subplot(2,2,2)
plot( config.k2/(2*pi)*config.c, 20*log10( abs( pw_b(:, 14) ) ) );
xlabel('f / Hz');
ylabel('magnitude / dB');
title('\alpha_{0} = \pi/2')

subplot(2,2,3)
plot( config.k2/(2*pi)*config.c, 20*log10( abs( pw_b(:, 28) ) ) );
xlabel('f / Hz');
ylabel('magnitude / dB');
title('\alpha_{0} = \pi')

subplot(2,2,4)
plot( config.k2/(2*pi)*config.c, 20*log10( abs( pw_b(:, 42) ) ) );
xlabel('f / Hz');
ylabel('magnitude / dB');
title('\alpha_{0} = 3\pi/2')


%plot phase response for chosen LS
figure;

subplot(2,2,1)
Q = unwrap( angle( pw_b(:,1) ) );
plot(config.k2/(2*pi)*config.c, Q);
xlabel('f / Hz');
ylabel('angle / rad');
title('LS bei \alpha_{0} = 0');

subplot(2,2,2)
Q = unwrap( angle( pw_b(:,14) ) );
plot(config.k2/(2*pi)*config.c, Q);
xlabel('f / Hz');
ylabel('angle / rad');
title('LS bei \alpha_{0} = \pi/2');

subplot(2,2,3)
Q = unwrap( angle( pw_b(:,28) ) );
plot(config.k2/(2*pi)*config.c, Q);
xlabel('f / Hz');
ylabel('angle / rad');
title('LS bei \alpha_{0} = \pi');

subplot(2,2,4)
Q = unwrap( angle( pw_b(:,42) ) );
plot(config.k2/(2*pi)*config.c, Q);
xlabel('f / Hz');
ylabel('angle / rad');
title('LS bei \alpha_{0} = 3\pi/2');