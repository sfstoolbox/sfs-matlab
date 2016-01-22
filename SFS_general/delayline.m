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
%                 conf.fracdelay.filter;
%                 if conf.fracdelay.filter~='zoh' 
%                   conf.fracdelay.order;                 
%                 conf.fracdelay.pre.method;
%                 if conf.fracdelay.pre.method=='resample'                   
%                   conf.fracdelay.pre.resample.method;
%                   conf.fracdelay.pre.resample.factor;
%                   if conf.fracdelay.pre.resample.method=='pm'
%                     conf.fracdelay.pre.resample.order;                   
%                 if conf.fracdelay.pre.method=='farrow'
%                 	conf.fracdelay.pre.farrow.Npol;
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
% Copyright (c) 2010-2016 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2016 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
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
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************

%% ===== Configuration ==================================================
fracdelay = conf.fracdelay;
Norder = conf.fracdelay.order;

%% ===== Computation =====================================================
% Check if the impulse response is given in SOFA conventions [M C N], or in
% usual [N C] convention, where
% M ... number of measurements
% C ... number of channels
% N ... number of samples
if ndims(sig)==3
  [M, C, samples] = size(sig);
  channels = M * C;
  % Reshape [M C N] => [N C*M], this will be redone at the end of the function
  sig = reshape(sig,[channels,samples])';
  reshaped = true;
else
  % Assume standard format [N C]
  [samples, channels] = size(sig);
  reshaped = false;
end

rfactor = 1.0;  % ratio of signal length and delayline length
switch fracdelay.pre.method
  case 'resample'
    % === Resample =====================================================
    % resample factor (1/stepsize of fractional delays)
    rfactor = fracdelay.pre.resample.factor;
    switch fracdelay.pre.resample.method
      case 'matlab'
        buffer = resample(sig,rfactor,1);
      case 'pm'
        % === Parks-McClellan linear phase FIR filter ===
        a = [1 1 0 0];
        f = [0.0 0.9/rfactor 1/rfactor 1.0];
        b = firpm(fracdelay.pre.resample.order,f,a);
        
        sig = reshape(sig, 1, channels*samples);
        sig = [sig; zeros(rfactor-1,channels*samples)];
        sig = reshape(sig, rfactor*samples, channels);
        
        buffer = filter(b,1,sig,[],1);
      otherwise
        disp('Delayline: Unknown resample method');
    end
  case 'farrow'
    % === Farrow-Structure ==============================================
    % Based on the assumption, that each coefficient h(n) of the fractional
    % delay filter can be expressed as a polynomial in d (frac. delay), i.e.
    %            _
    %           \  NPol
    % h_d(n) ~=  >      c_m(n) d^m
    %           /_ m=0
    %
    % For some Filter design methods, e.g. Lagrange Interpolators, this
    % perfectly possible. For other, a uniform grid of test delays d_q is
    % used to fit the polynomials to the desired coefficient(n) find a set
    % polynomial which approximates each coefficient of the desired filter.
    % This structure allows to perform the convolution independently from the
    % delay and reuse the results of the filter for different delays.
    %                          _
    %                          \  NPol
    % y(n) = h_d(n) * x(n) ~=   >      ( c_m(n)*x(n) ) d^m
    %                          /_ m=0
    %
    % The above representation shows that the convolution of the input signal x
    % can be performed by first convolving c_m and x and incorporating the delay
    % d afterwards.

    % number of filter coefficients (length of filter)
    % Norder = fracdelay.order;  
    % number of parallel filters, i.e. order of polynomial + 1
    % Nfilter = fracdelay.pre.farrow.Npol+1;
    
    to_be_implemented(mfilename);
  case 'none'
    buffer = sig;
  otherwise
    disp('%s: \"%s\" is an unknown pre-processing method for delay line', ...
      fracdelay.pre.method, upper(mfilename));
end

%% ===== Fractional Delay ================================================

dt = rfactor.*dt;  % resampled delays 
samples = rfactor.*samples;  % length of resampled signals

if strcmp( fracdelay.pre.method, 'farrow')  
  to_be_implemented(mfilename);
else  % There is no post processing stage if the Farrow Structure used
  % === Post Processing ====================================================
  a = ones(1, channels);  % denominator of fractional delay filter
  switch fracdelay.filter
    case 'zoh'
      % === Zero-Order-Hold (Integer Delays) ===============================
      idt = ceil(dt);  % round up to next integer delay
      b = ones(1, channels);
    case 'lagrange'
      % ==== Lagrange Polynomial Interpolator ==============================
      if mod(Norder,2) == 0
        idt = round(dt);  % round delay for even order
      else
        idt = floor(dt);  % floor delay for odd order
      end
      fdt = dt - idt;  % fractional part of delays
      b = lagrange_filter(Norder, fdt);
    case 'thiran'
      % ==== Thiran's Allpass Filter for Maximally Flat Group Delay ========
      idt = round(dt);  % integer part of delays
      fdt = dt - idt;  % fractional part of delays
      [b, a] = thiran_filter(Norder, fdt);
    case 'least_squares'
      % ==== Least Squares Interpolation Filter ============================
      idt = floor(dt);  % integer part of delays
      fdt = dt - idt;  % fractional part of delays
      b = zeros(Norder+1, channels);
      for cdx=1:channels
        b(:,cdx) = general_least_squares(Norder+1,fdt(cdx),0.90);
      end
    otherwise
      error('%s: \"%s\" is an unknown fractional delay filter', ...
        upper(mfilename), fracdelay.filter);
  end
  
  for cdx=1:channels
    buffer(:,cdx) = filter(b(:,cdx),a(:,cdx),buffer(:,cdx));
  end  
end

%% ===== Integer Delay ================================================
% Handling of too long delay values (returns vector of zeros)
idt(abs(idt) > samples) = samples;

% Handle positive or negative delays
for cdx=1:channels
  if idt(cdx)>=0
    buffer(:,cdx) = [zeros(idt(cdx),1); weight(cdx)*buffer(1:end-idt(cdx),cdx)];
  else
    buffer(:,cdx) = [weight(cdx)*buffer(-idt(cdx)+1:end,cdx); zeros(-idt(cdx),1)];
  end
end
sig = buffer(1:rfactor:samples,:);

%%
% Undo reshaping [N M*C] => [M C N]
if reshaped
  sig = reshape(sig',[M C size(sig,1)]);
end