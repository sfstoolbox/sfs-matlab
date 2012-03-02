function [s] =  delayline(s,dt,weight,conf)
%DELAYLINE implements a (fractional) weighting delay line
%
%   Usage: s = delayline(s,dt,weight,conf)
%          
%
%
%   Input parameters:
%       s       - input signal (vector)
%       dt      - delay (samples)
%       weight  - weighting factor
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output:
%       s    - delayed signal

% AUTHOR: Sascha Spors
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$



%% ===== Checking of input  parameters ==================================


%% ===== Configuration ==================================================
fracdelay = conf.usefracdelay;
method = conf.fracdelay_method;

rfactor=10;         % resample factor
Lls=30;             % length of least-squares factional delay filter

%% ===== Computation =====================================================    
if(fracdelay)
    
    conf2.fracdelay=0;
    conf2.fracdelay_method='';

    switch method

    case 'resample'
       s2 = resample(s,rfactor,1);
       s2 = delayline(s2,rfactor*dt,weight,conf2);
       s = resample(s2,1,rfactor);

    case 'least_squares'
        idt=floor(dt);
        s = delayline(s,idt,weight,conf2);

        if(abs(dt-idt)>0)
            h = hgls2(Lls,dt-idt,0.90);
            s = conv(s,h);
            s = s(Lls/2:end-Lls/2);
        end

    otherwise
        disp('Delayline: Unknown fractional delay method');

    end
    
else
% from here on integer delays are considered
    idt=round(dt);

    if(idt==0)
        s = weight*s;
        return;
    end

    if(idt > 0)
        s = [zeros(1,idt) weight*s(1:end-idt)];
    else
        s = [weight*s(-idt+1:end) zeros(1,-idt)];
    end

end

