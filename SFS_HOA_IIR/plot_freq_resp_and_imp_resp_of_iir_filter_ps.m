%plot magnitude response for chosen LS
figure;

subplot(2,2,1)
hold on
[y1]=freqz(OUT(:,1),1,config.f2,config.fs);
plot(config.f2,db(y1))
plot( config.k2/(2*pi)*config.c, 20*log10( abs( E(:, 1) ) ),'r:' );
xlabel('frequency');
ylabel('magnitude / dB');
title('\alpha_{0} = 0')
grid on

subplot(2,2,2)
hold on
[y2]=freqz(OUT(:,14),1,config.f2,config.fs);
plot(config.f2,db(y2))
plot( config.k2/(2*pi)*config.c, 20*log10( abs( E(:, 14) ) ),'r:' );
xlabel('frequency');
ylabel('magnitude / dB');
title('\alpha_{0} = \pi/2')
grid on

subplot(2,2,3)
hold on
[y3]=freqz(OUT(:,28),1,config.f2,config.fs);
plot(config.f2,db(y3))
plot( config.k2/(2*pi)*config.c, 20*log10( abs( E(:, 28) ) ),'r:' );
xlabel('frequency');
ylabel('magnitude / dB');
title('\alpha_{0} = \pi')
grid on

subplot(2,2,4)
hold on
[y4]=freqz(OUT(:,42),1,config.f2,config.fs);
plot(config.f2,db(y4))
plot( config.k2/(2*pi)*config.c, 20*log10( abs( E(:, 42) ) ),'r:' );
xlabel('frequency');
ylabel('magnitude / dB');
title('\alpha_{0} = 3\pi/2')
grid on

%plot phase response for chosen LS
figure;

subplot(2,2,1)
hold on
P = unwrap( -angle( y1(:,1) ) );
plot(config.f2, P);
Q = unwrap( angle( E(:,1) ) );
plot(config.f2, Q,'r:');
xlabel('f / Hz');
ylabel('angle / rad');
title('LS bei \alpha_{0} = 0');
grid on

subplot(2,2,2)
hold on
P = unwrap( -angle( y2(:,1)*180/pi ) );
plot(config.f2, P);
Q = unwrap( angle( E(:,14) ) );
plot(config.f2, Q,'r:');
xlabel('f / Hz');
ylabel('angle / rad');
title('LS bei \alpha_{0} = \pi/2');
grid on

subplot(2,2,3)
hold on
P = unwrap( -angle( y3(:,1)*180/pi ) );
plot(config.f2, P);
Q = unwrap( angle( E(:,28) ) );
plot(config.f2, Q,'r:');
xlabel('f / Hz');
ylabel('angle / rad');
title('LS bei \alpha_{0} = \pi');
grid on

subplot(2,2,4)
hold on
P = unwrap( -angle( y4(:,1)*180/pi ) );
plot(config.f2, P);
Q = unwrap( angle( E(:,42) ) );
plot(config.f2, Q,'r:');
xlabel('f / Hz');
ylabel('angle / rad');
title('LS bei \alpha_{0} = 3\pi/2');
grid on



%plot impulse response
figure;

subplot(2,2,1)
hold on
plot(OUT(:,1));
xlabel('samples');
ylabel('amplitude');
title('Impulse response, LS bei \alpha_{0} = 0');


subplot(2,2,2)
hold on
plot(OUT(:,14));
xlabel('samples');
ylabel('amplitude');
title('Impulse response, LS bei \alpha_{0} = \pi/2');


subplot(2,2,3)
hold on
plot(OUT(:,28));
xlabel('samples');
ylabel('amplitude');
title('Impulse response, LS bei \alpha_{0} = \pi');


subplot(2,2,4)
hold on
plot(OUT(:,42));
xlabel('samples');
ylabel('amplitude');
title('Impulse response, LS bei \alpha_{0} = 3\pi/2');
