function sig = delayline(sig,dt,weight,conf)
%DELAYLINE implements a (fractional) delay line with weights
%
%   Usage: sig = delayline(sig,dt,weight,conf)
%
%   Input parameter:
%       sig     - input signal (vector), can be in the form of [N C], or
%                 [M C N], where
%                     N ... samples
%                     C ... channels (most probably 2)
%                     M ... number of measurements
%                 If the input is [M C N], the length of dt and weight has to be
%                 1 or M*C. In the last case the first M entries in dt are
%                 applied to the first channel and so on.
%       dt      - delay / samples
%       weight  - amplitude weighting factor
%       conf    - configuration struct (see SFS_config).
%                 Used settings are:
%                     conf.usefracdelay;
%                     conf.fracdelay_method; (only if conf.usefracdelay==true)
%
%   Output parameter:
%       sig     - delayed signal
%
%   DELAYLINE(sig,dt,weight,conf) implementes a delayline, that delays the given
%   signal by dt samples and applies an amplitude weighting factor. The delay is
%   implemented as integer delay or fractional delay filter, see description of
%   conf input parameter.
%
%   See also: get_ir, driving_function_imp_wfs

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Team                                   *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Configuration ==================================================
usefracdelay = conf.usefracdelay;


%% ===== Computation =====================================================
% Check if the impulse response is given in SOFA conventions [M C N], or in
% usual [N C] convention, where
% M ... number of measurements
% C ... number of channels
% N ... number of samples
if ndims(sig)==3
    [M C samples] = size(sig);
    channels = M * C;
    % Reshape [M C N] => [N C*M], this will be redone at the end of the function
    sig = reshape(sig,[channels,samples])';
    reshaped = true;
else
    % Assume standard format [N C]
    [samples channels] = size(sig);
    reshaped = false;
end
% If only single valued time delay and weight is given, create vectors
if channels>1 && length(dt)==1, dt=repmat(dt,[1 channels]); end
if channels>1 && length(weight)==1, weight=repmat(weight,[1 channels]); end

if usefracdelay

    % Additional configuration
    fracdelay_method = conf.fracdelay_method;
    rfactor = 100; % resample factor (1/stepsize of fractional delays)
    Lls = 30;      % length of least-squares factional delay filter

    % Defining a temporary conf struct for recursive calling of delayline
    conf2.usefracdelay = false;

    switch fracdelay_method
    case 'resample'
       sig2 = resample(sig,rfactor,1);
       sig2 = delayline(sig2,rfactor.*dt,weight,conf2);
       sig = resample(sig2,1,rfactor);

    case 'least_squares'
        idt = floor(dt);
        sig = delayline(sig,idt,weight,conf2);
        if abs(dt-idt)>0
            for ii=1:channels
                h = general_least_squares(Lls,dt(ii)-idt(ii),0.90);
                tmp = convolution(sig(:,ii),h);
                sig(:,ii) = tmp(Lls/2:end-Lls/2);
            end
        end

    case 'interp1'
        idt = floor(dt);
        sig = delayline(sig,idt,weight,conf2);
        t = (1:samples)';
        for ii=1:channels
            sig(:,ii) = interp1(t,sig(:,ii),-(dt(ii)-idt(ii))+t,'spline');
        end

    otherwise
        disp('Delayline: Unknown fractional delay method');
    end

else
    % From here on integer delays are considered
    idt = round(dt);

    % Handling of too long delay values (returns vector of zeros)
    idt(abs(idt)>samples) = samples;

    % Handle positive or negative delays
    for ii=1:channels
        if idt(ii)>=0
            sig(:,ii) = [zeros(idt(ii),1); weight(ii)*sig(1:end-idt(ii),ii)];
        else
            sig(:,ii) = [weight(ii)*sig(-idt(ii)+1:end,ii); zeros(-idt(ii),1)];
        end
    end
end

% Undo reshaping [N M*C] => [M C N]
if reshaped
    sig = reshape(sig',[M C samples]);
end
