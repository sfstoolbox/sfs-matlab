function sig = delayline(sig,dt,weight,conf)
%DELAYLINE implements a (fractional) weighting delay line
%   Usage: sig = delayline(sig,dt,weight,conf)
%          sig = delayline(sig,dt,weight)
%          
%   Input parameter:
%       sig     - input signal (vector)
%       dt      - delay (samples)
%       weight  - weighting factor
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameter:
%       sig     - delayed signal

% AUTHOR: Sascha Spors
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 4;
error(nargchk(nargmin,nargmax,nargin))
if nargin==nargmax-1
    conf = SFS_config;
end


%% ===== Configuration ==================================================
usefracdelay = conf.usefracdelay;
fracdelay_method = conf.fracdelay_method;

rfactor=100; % resample factor (1/stepsize of fractional delays)
Lls=30;      % length of least-squares factional delay filter


%% ===== Computation =====================================================    
if(usefracdelay)
 
    % Defining a temporary conf struct for recursive calling of delayline
    conf2.usefracdelay=0;
    conf2.fracdelay_method='';

    switch fracdelay_method
    case 'resample'
       sig2 = resample(sig,rfactor,1);
       sig2 = delayline(sig2,rfactor*dt,weight,conf2);
       sig = resample(sig2,1,rfactor);

    case 'least_squares'
        idt=floor(dt);
        sig = delayline(sig,idt,weight,conf2);
        if(abs(dt-idt)>0)
            % FIXME: add this function to Toolbox, or display URL for download
            h = hgls2(Lls,dt-idt,0.90);
            %[IP,wprot] = iniheq2(Lls,0.90);
            %h = heqrip2(Lls,dt-idt,wprot,IP);
            sig = conv(sig,h);
            sig = sig(Lls/2:end-Lls/2);
        end

    case 'interp1'
        idt=floor(dt);
        sig = delayline(sig,idt,weight,conf2);
        t=1:length(sig);
        sig = interp1(t,sig,-(dt-idt)+t,'spline');
        
    otherwise
        disp('Delayline: Unknown fractional delay method');
    end
    
else
    % from here on integer delays are considered
    idt=round(dt);
    % handle positive or negative delays
    if idt>=0
        sig = [zeros(1,idt) weight*sig(1:end-idt)];
    else
        sig = [weight*sig(-idt+1:end) zeros(1,-idt)];
    end
end
