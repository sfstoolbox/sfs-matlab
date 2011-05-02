function [OUT, X, F, Gd, F2] = HOA_filter_firintpolc(E, config)
    
OUT = zeros(config.N + 1, config.nLS);

for n = 0:1:config.nLS-1

    disp([n config.nLS-1]); % display computation time
    [h0] = firintpolc(config.N, config.f3, E(:,n+1), config.g); %calculation of impulse responses
    OUT(:,n+1) =  h0; % write impulse responses in matrix OUT
    
end

[X,F] = freqz( OUT(:,config.cLS), 1, config.n, 40000); % calculation of frequency response
[Gd,F2] = grpdelay( OUT(:,config.cLS), 1, config.n, config.fs); % group delay of FIR Filter


% Fourier Transformation of impulse response to get frequency response

% L = config.N;                  % Length of signal
% NFFT = 2^nextpow2(L); % Next power of 2 from length of L
% Y = fft(OUT(:,config.cLS),NFFT);
% f = config.fs/2*linspace(0,1,NFFT/2+1);
% 
% % Plot single-sided amplitude spectrum.
% figure
% hold on
% plot(f,db(abs(Y(1:NFFT/2+1))))
% plot(config.f2, db( abs( E(:,config.cLS) ) ), 'r')
% hold off


