function delay = delay_line_fir_pw(config)

% x and y coordinate of wave vector
kx = cos(config.al_pw);
ky = sin(config.al_pw);

% calculation of delay for every LS
delay = zeros(1,length(config.x0));

    for n=0:length(config.x0)-1 %loop over all loudspeakers
    
        delay(n+1) = 1/config.c*(kx*config.x0(1,n+1)+ky*config.x0(2,n+1));
    
    end

%delay in time domain 
delay = delay + max(abs(delay)); 

% delay in samples
delay = -round(delay.*config.fs); 

end