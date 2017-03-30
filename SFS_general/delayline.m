function [sig,delay_offset] = delayline(sig,dt,weight,conf)
%DELAYLINE implements a (fractional) delay line with weights
%
%   Usage: [sig,delay_offset] = delayline(sig,dt,weight,conf)
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
%       dt      - delay / s
%       weight  - amplitude weighting factor
%       conf    - configuration struct (see SFS_config).
%
%   Output parameter:
%       sig             - delayed signal
%       delay_offset    - additional delay / s
%                         This is added by the fractional delayline filters to
%                         all channels. For integer delays this is 0.
%
%   DELAYLINE(sig,dt,weight,conf) implementes a delayline, that delays the given
%   signal by dt samples and applies an amplitude weighting factor. The delay is
%   implemented as integer delay or fractional delay filter, see delayline
%   section in SFS_config for possible settings. As default setting an integer
%   delayline is used.
%
%   See also: get_ir, driving_function_imp_wfs

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
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


%% ===== Configuration ===================================================
% Check for old configuration
if isfield(conf, 'usefracdelay')
    error(['%s: conf.usefracdelay is deprecated, please use conf.delayline', ...
           ' instead. See SFS_config for details.'],upper(mfilename));
end
fs = conf.fs;
delay = conf.delayline;


%% ===== Preparation =====================================================
% --- Reshape signals ---
% Check if the signal is an impulse response given in SOFA conventions [M C N],
% or in usual [N C] convention, where
% M ... number of measurements
% C ... number of channels
% N ... number of samples
if ndims(sig)==3
    [M,C,samples] = size(sig);
    channels = M * C;
    % Reshape [M C N] => [N C*M], this will be redone at the end of the function
    sig = reshape(sig,[channels,samples])';
    reshaped = true;
else
    % Assume standard format [N C]
    [samples,channels] = size(sig);
    reshaped = false;
end


%% ===== Resampling ======================================================
% The resampling is applied independently from the actual fractional/integer
% delay handling performed in the next step. The resampling is redone at the end
% of the file.
% If resampling is used together with the integer delay filter, this is already
% a usage of fractional delay due to the upsampling.
%
switch delay.resampling
case 'none'
    rfactor = 1.0;
    delay_offset = 0;
case 'matlab'
    rfactor = delay.resamplingfactor;
    delay_offset = 0;
    sig = resample(sig,rfactor,1);
case 'pm'
    % === Parks-McClellan linear phase FIR filter ===
    rfactor = delay.resamplingfactor;
    rfilt = pm_filter(rfactor*delay.resamplingorder, 0.9/rfactor, ...
      1/rfactor);
    delay_offset = delay.resamplingorder*rfactor/2;

    sig = reshape(sig,1,channels*samples);
    sig = [sig; zeros(rfactor-1,channels*samples)];
    sig = reshape(sig,rfactor*samples,channels);

    sig = convolution(rfactor*rfilt, sig);
otherwise
    error('%s: "%s": unknown resampling method',upper(mfilename), ...
        delay.resampling);
end


%% ===== Expansion of signals, delays or weights =========================
% --- Expand channels
if channels==1
    channels = max(length(dt),length(weight));
    sig = repmat(sig,[1 channels]);
end

% --- Expand dt and weight ---
% If only single valued time delay and weight is given, create vectors
if channels>1 && length(dt)==1, dt=repmat(dt,[1 channels]); end
if channels>1 && length(weight)==1, weight=repmat(weight,[1 channels]); end


%% ===== Conversion to integer delay =====================================
dt = dt.*rfactor.*fs;  % resampled delay / samples
samples = rfactor.*samples;  % length of resampled signals
switch delay.filter
case 'integer'
    % === Integer delays ===
    idt = round(dt);  % round to nearest integer delay
    delay_offset = delay_offset + 0;
case 'zoh'
    % === Zero-order hold ===
    idt = ceil(dt);  % round to next larger integer delay
    delay_offset = delay_offset + 0;
case 'lagrange'
    % === Lagrange polynomial interpolator ===
    if iseven(delay.filterorder)
        idt = round(dt);  % round delay for even order
    else
        idt = floor(dt);  % floor delay for odd order
    end
    fdt = dt - idt;  % fractional part of delays
    b = lagrange_filter(delay.filterorder,fdt);
    a = ones(1,channels);
    delay_offset = delay_offset + floor(delay.filterorder/2);
case 'thiran'
    % === Thiran's allpass filter for maximally flat group delay ===
    idt = round(dt);  % integer part of delays
    fdt = dt - idt;  % fractional part of delays
    [b,a] = thiran_filter(delay.filterorder,fdt);
    delay_offset = delay_offset + delay.filterorder;
case 'least_squares'
    % ==== Least squares interpolation filter ===
    idt = floor(dt);  % integer part of delays
    fdt = dt - idt;  % fractional part of delays
    b = zeros(delay.filterorder+1,channels);
    for ii=1:channels
        b(:,ii) = general_least_squares(delay.filterorder+1,fdt(ii),0.90);
    end
    a = ones(1,channels);
    delay_offset = delay_offset + floor(delay.filterorder/2);
case 'farrow'
    % === Farrow-structure ===
    % Based on the assumption, that each coefficient h(n) of the fractional
    % delay filter can be expressed as a polynomial in d (frac. delay), i.e.
    %            __
    %           \  NPol
    % h_d(n) ~=  >      c_m(n) d^m
    %           /__m=0
    %
    % For some Filter design methods, e.g. Lagrange Interpolators, this is
    % perfectly possible. For other, a uniform grid of test delays d_q is
    % used to fit the polynomials to the desired coefficient(n) find a set
    % polynomial which approximates each coefficient of the desired filter.
    % This structure allows to perform the convolution independently from
    % the delay and reuse the results of the filter for different delays.
    %                           __
    %                          \  NPol
    % y(n) = h_d(n) * x(n) ~=   >      ( c_m(n)*x(n) ) d^m
    %                          /__m=0
    %
    % The above representation shows that the convolution of the input
    % signal x can be performed by first convolving c_m and x and
    % incorporating the delay d afterwards.
    %
    % number of parallel filters, i.e. order of polynomial + 1
    % Nfilter = delay.filternumber;
    to_be_implemented(mfilename);
otherwise
    error('%s: \"%s\" is an unknown delayline filter', ...
        upper(mfilename),delay.filter);
end
% Apply filter if needed
if exist('a','var') && exist('b','var')
    for ii=1:channels
        sig(:,ii) = filter(b(:,ii),a(:,ii),sig(:,ii));
    end
end


%% ===== Integer delayline ===============================================
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


%% ===== Postprocessing ==================================================
% --- Downsampling ---
if rfactor~=1
    sig = sig(1:rfactor:samples,:);
    delay_offset = delay_offset ./ rfactor;
end
% --- Undo reshape ---
% [N M*C] => [M C N]
if reshaped
    % C might have changed due to replication of single-channel input
    sig = reshape(sig',M,[],size(sig,1));
end
% --- delay_offset in seconds ---
delay_offset = delay_offset / fs;
