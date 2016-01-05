function delayline = delayline_write(sig,conf)
%DELAYLINE implements a (fractional) delay line with weights
%
%   Usage: sig = delayline_write(sig,[conf])
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
%       conf    - configuration struct (see SFS_config).
%                 Used settings are:
%                     conf.usefracdelay;
%                     conf.fracdelay_method; (only if conf.usefracdelay==true)
%
%   Output parameter:
%       sig     - data stream representing the delayline
%
%   DELAYLINE(sig,dt,weight,conf) implementes a delayline, that delays the given
%   signal by dt samples and applies an amplitude weighting factor. The delay is
%   implemented as integer delay or fractional delay filter, see description of
%   conf input parameter.
%
%   See also: get_ir, driving_function_imp_wfs

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
    Nh = fracdelay.length;  % number of filter coefficients (length of filter)    
    switch fracdelay.filter
      case 'lagrange'
        % ==== Lagrange Polynomial Interpolator ==============================
        NFilter = Nh;  % number of parallel filters
        % each row is a polynom in d, each column is a filter
        c = lagrange_polynoms(0:(Nh-1));
      otherwise
        disp('Delayline: Filter not implemented in farrow structure');
    end
    delayline = zeros(NFilter*samples, channels);
    % store output of each parallel filter as interleaved data
    for ndx=1:NFilter
      delayline(ndx:NFilter:end, :) = filter(c(:,ndx), 1, sig, [], 1);
    end
    rfactor = NFilter;
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
