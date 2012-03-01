function [s] =  delayline(s,dt,weight,conf)
%DELAYLINE implements a (fractional) weighting delay line
%
%   Usage: s = delayline(s,dt,weight,hf)
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
fracdelay = conf.fracdelay;
method = conf.fracdelay_method;

%% ===== Computation =====================================================    
if(fracdelay==0)
    idt=round(dt);

    if(idt==0)
        return;
    end

    if(idt > 0)
        s = [zeros(1,idt) weight*s(1:end-idt)];
    else
        s = [weight*s(-idt+1:end) zeros(1,-idt)];
    end
return;    
end

% from here on fractional delays are considered
conf2.fracdelay=0;
conf2.fracdelay_method='';

switch method
   
    case 'resample'
       rfactor=10;
       s2 = resample(s,rfactor,1);
       s2 = delayline(s2,rfactor*dt,weight,conf2);
       s = resample(s2,1,rfactor);
       
    case 'lagrange'
        idt=floor(dt);
        s = delayline(s,idt,weight,conf2);
        
        if(abs(dt-idt)>0)
            L=30;
            %h = hlagr2(L,dt-idt);
            h = hgls2(L,dt-idt,0.90);
            s = conv(s,h);
            s = s(L/2:end-L/2);
            
            %s = [s(L/2:end) zeros(1,L/2)];
        end
    
    otherwise
        disp('Delayline: Unknown fractional delay method');
    
end

