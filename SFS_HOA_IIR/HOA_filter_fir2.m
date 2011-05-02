% OUT holds n+1 coefficients of the digital FIR Filter for n LS

function [OUT, X, F, Gd, F2] = HOA_filter_fir2(E, config, delay)

% calculation of FIR Filter coefficients for all LS 
OUT = zeros(config.N + 1, config.nLS);

for n = 0:1:config.nLS-1
     
    % filter for every loudspeaker
    b = fir2(config.N, config.f3, abs(E(:,n+1).') ); 
    
    % impulse responses
    OUT(:,n+1) =  b;

    % shift impulse response to get the desired delay 
    OUT(:,n+1) = circshift(OUT(:,n+1),[(delay(n+1)) 0]); 
    
end

% frequency response of FIR Filter for choosen LS
[X,F] = freqz( OUT(:,config.cLS), 1, config.n, config.fs); 

% group delay of FIR Filter
[Gd,F2] = grpdelay( OUT(:,config.cLS), 1, config.n, config.fs);






