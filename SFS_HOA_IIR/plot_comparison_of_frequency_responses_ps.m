figure;

subplot(2,2,1)
plot( config.k2/(2*pi)*config.c, 20*log10( abs(ps_b(:, 1) ) ) );
hold on
plot( config.k2/(2*pi)*config.c, 20*log10( abs(ps_b_r(:, 1) ) ), 'r:');
xlabel('f / Hz');
ylabel('magnitude / dB');
title('LS 1');

subplot(2,2,2)
plot( config.k2/(2*pi)*config.c, 20*log10( abs(ps_b(:, 14) ) ) );
hold on
plot( config.k2/(2*pi)*config.c, 20*log10( abs(ps_b_r(:, 14) ) ),'r:' );
xlabel('f / Hz');
ylabel('magnitude / dB');
title('LS 14');

subplot(2,2,3)
plot( config.k2/(2*pi)*config.c, 20*log10( abs(ps_b(:, 28) ) ) );
hold on
plot( config.k2/(2*pi)*config.c, 20*log10( abs(ps_b_r(:, 28) ) ),'r:');
xlabel('f / Hz');
ylabel('magnitude / dB');
title('LS 28');

subplot(2,2,4)
plot( config.k2/(2*pi)*config.c, 20*log10( abs(ps_b(:, 42) ) ) );
hold on
plot( config.k2/(2*pi)*config.c, 20*log10( abs(ps_b_r(:, 42) ) ),'r:' );
xlabel('f / Hz');
ylabel('magnitude / dB');
title('LS 42');