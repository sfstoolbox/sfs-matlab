function delayline = delayline_write(sig,conf)
%DELAYLINE_WRITE writes data into a (fractional) delay line
%
%   Usage: sig = delayline_write(sig,conf)
%
%   Input parameter:
%       sig     - input signal (vector), can be in the form of [N C], or
%                 [M C N], where
%                     N ... samples
%                     C ... channels (most probably 2)
%                     M ... number of measurements
%       conf    - configuration struct (see SFS_config).
%                 Used settings are:
%                   conf.fracdelay.pre.method;
%                   
%                   (only if conf.fracdelay.pre.method=='resample')
%                   conf.fracdelay.pre.resample.method
%                   conf.fracdelay.pre.resample.factor
%                   
%                   (only if conf.fracdelay.pre.method=='farrow')
%                   conf.fracdelay.order;
%                   fracdelay.pre.farrow.Npol
%
%   Output parameter:
%       delayline - data structure representing the delayline
%
%   DELAYLINE_WRITE(sig,conf) implements the pre-processing stage of a 
%   fractional delay line for a given signal. This processing stage generates
%   a data structure which is indepedent of the actual delays applied to the 
%   signal. The data structure may than be used as an input to delayline_read
%   in order to extract multiple delayed versions of the signal.
%
%   See also: get_ir, driving_function_imp_wfs, delayline_read

%*****************************************************************************
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
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
        delayline = resample(sig,rfactor,1);
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
    Norder = fracdelay.order;  
    % number of parallel filters, i.e. order of polynomial + 1
    Nfilter = fracdelay.pre.farrow.Npol+1;
    
    if strcmp(fracdelay.filter, 'lagrange')
      % ==== Lagrange Polynomial Interpolator ===============================
      if Nfilter ~= Norder+1
        error( ['%s: order of lagrange filter has be to equal to the number', ...
          ' of parallel filters of farrow structure'] ...
          , upper(mfilename));
      end
      % each row is a polynom in d, each column is a filter
      c = lagrange_polynomials(((1-Norder)/2:(1+Norder)/2));      
    else     
      d = (0:(Nfilter-1))./Nfilter;  % uniform grid of delays to fit polynomials
      d(1) = d(1)+1E-5; % prevent some issues with d=0.0;
      h = zeros(Norder+1, Nfilter); % prototype filters;
      switch fracdelay.filter
       case 'least_squares'
        % ==== General Least Squares Method =================================
        for hdx=1:Nfilter
          h(:,hdx) = general_least_squares(Norder+1,d(hdx),0.90);
        end
      otherwise
        disp('Delayline: Filter not implemented in farrow structure');
      end
      % fit polynomials using least squares approximation
      c = zeros(Norder+1,Nfilter);
      for ndx=1:Norder+1
        c(ndx,:) = polyfit(d,h(ndx,:),Nfilter-1);
      end      
    end
    delayline = zeros(Nfilter*samples, channels);
    % store output of each parallel filter as interleaved data
    for ndx=1:Nfilter
      delayline(ndx:Nfilter:end, :) = filter(c(:,ndx), 1, sig, [], 1);
    end
    rfactor = Nfilter;
  case 'none'
    delayline = sig;
  otherwise
    disp('Delayline: Unknown Pre-Processing method for delay line');
end

%%
% Undo reshaping [N M*C] => [M C N]
if reshaped
  delayline = reshape(delayline',[M C rfactor*samples]);
end
