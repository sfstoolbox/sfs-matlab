function sig = delayline(sig,dt,weight,conf)
%DELAYLINE implements a (fractional) weighting delay line
%
%   Usage: sig = delayline(sig,dt,weight,[conf])
%
%   Input parameter:
%       sig     - input signal (vector)
%       dt      - delay (samples)
%       weight  - amplitude weighting factor
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameter:
%       sig     - delayed signal
%
%   DELAYLINE(sig,dt,weight,conf) implementes a delayline, that delays the given
%   signal by dt samples and applies a amplitude weighting factor.
%   The delay is implemented as integer delays or a fractional delay
%   filter, see SFS_config.
%
%   see also: driving_function_imp_wfs_25d

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax)
if nargin==nargmax-1
    conf = SFS_config;
end


%% ===== Configuration ==================================================
usefracdelay = conf.usefracdelay;
fracdelay_method = conf.fracdelay_method;

rfactor = 100; % resample factor (1/stepsize of fractional delays)
Lls = 30;      % length of least-squares factional delay filter


%% ===== Computation =====================================================
if(usefracdelay)

    % Defining a temporary conf struct for recursive calling of delayline
    conf2.usefracdelay = 0;
    conf2.fracdelay_method = '';

    switch fracdelay_method
    case 'resample'
       sig2 = resample(sig,rfactor,1);
       sig2 = delayline(sig2,rfactor*dt,weight,conf2);
       sig = resample(sig2,1,rfactor);

    case 'least_squares'
        if ~exist('hgls2')
            error(['%s: the least_squares methods needs the hgls2 function ',...
                'which you have to look for in the web ;)']);
        end
        idt = floor(dt);
        sig = delayline(sig,idt,weight,conf2);
        if(abs(dt-idt)>0)
            h = hgls2(Lls,dt-idt,0.90);
            sig = conv(sig,h);
            sig = sig(Lls/2:end-Lls/2);
        end

    case 'interp1'
        idt = floor(dt);
        sig = delayline(sig,idt,weight,conf2);
        t = 1:length(sig);
        sig = interp1(t,sig,-(dt-idt)+t,'spline');

    otherwise
        disp('Delayline: Unknown fractional delay method');
    end

else
    % from here on integer delays are considered
    idt = round(dt);
    
    % handling of too long delay values (returns vector of zeros)
    if(abs(idt)>length(sig))
        idt=length(sig);
    end
    
    % handle positive or negative delays
    if idt>=0
        sig = [zeros(1,idt) weight*sig(1:end-idt)];
    else
        sig = [weight*sig(-idt+1:end) zeros(1,-idt)];
    end
end
