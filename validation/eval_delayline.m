% Evaluates the delayline implementation

% AUTHOR: Sascha Spors
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$

clear all;

%% ===== Configuration ==================================================

% delays to evaluate
%dt=[-1 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
%dt=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
%dt=[-5 -2.5 0 2.5 5];
dt=1*[0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1];

% parameters for delayline
conf.fracdelay = 1;
conf.fracdelay_method = 'least_squares';

% length of input signal
L=256;

% create frequency axis
w = (0:1:(L-1))/L; 
wpi = w*pi;
wpi2=wpi(2:L);

% set up input signal
insig = zeros(1,L);
insig(L/2) = 1;


%% ===== Computation =====================================================
for n=1:length(dt)
    outsig(:,n) = delayline(insig,dt(n),1,conf);

    H(:,n) = freqz(outsig(:,n),1,wpi);
    magresp(:,n) = abs(H(:,n));
    uwphase(:,n)=-unwrap(angle(H(:,n)));
    
    phasdel(:,n) = uwphase(2:L,n)./wpi2';
end



%% ===== Plotting =====================================================
% setup legend and axis
t=1:L;
t=t-L/2;

% phase delay
figure;
plot(wpi2/pi,phasdel-(L/2)+1);
ylabel('phase delay');
xlabel('normalized frequency');
grid on;

% magnitude response
figure;
plot(wpi/pi,magresp);
ylabel('magnitude');
xlabel('normalized frequency');
grid on;

% impluse response
figure;
plot(t(L/2-10:L/2+10),outsig(L/2-10:L/2+10,:));
ylabel('amplitude');
xlabel('samples');
grid on;
