function sig = delayline_read(delayline,dt,weight,conf)
%DELAYLINE implements a (fractional) delay line with weights
%
%   Usage: sig = delayline_read(sig,dt,weight,[conf])
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
Nh = fracdelay.length;  % length fractional delay filter

%% ===== Computation =====================================================
% Check if the impulse response is given in SOFA conventions [M C N], or in
% usual [N C] convention, where
% M ... number of measurements
% C ... number of channels
% N ... number of samples
if ndims(delayline)==3
  [M, C, samples] = size(delayline);
  channels = M * C;
  % Reshape [M C N] => [N C*M], this will be redone at the end of the function
  delayline = reshape(delayline,[channels,samples])';
  reshaped = true;
else
  % Assume standard format [N C]
  [samples, channels] = size(delayline);
  reshaped = false;
end
% If only single valued time delay and weight is given, create vectors
if channels>1 && length(dt)==1, dt=repmat(dt,[1 channels]); end
if channels>1 && length(weight)==1, weight=repmat(weight,[1 channels]); end

%% ===== Fractional Delay ================================================

rfactor = 1.0;  % ratio of signal length and delayline length
idt = floor(dt);  % integer part of delays
fdt = dt - idt;  % fractional part of delays

% Stuff depending on the pre processing method
switch fracdelay.pre.method
  case 'resample'
    % === Resample =====================================================
    % Assuming a resampled delay line
    rfactor = fracdelay.pre.resample.factor;
    dt = rfactor.*dt;
    idt = floor(dt);
    fdt = dt - idt;
  case 'farrow'
    % === Farrow-Structre ==============================================
    % number of parallel filters, i.e. order of polynomial - 1
    Nfilter = fracdelay.pre.farrow.Npol+1;
    samples = samples/Nfilter;
    tmp = zeros(samples, channels);
    for cdx=1:channels
	  % shorter, faster way of matlab's polyval
      d_vec = fdt(cdx).^(Nfilter-1:-1:0);
      tmp(:, cdx) = ...
        d_vec*reshape( delayline(:, cdx), [Nfilter, samples] );
    end
    delayline = tmp;
  case 'none'
  otherwise
    disp('Delayline: Unknown Pre-Processing method for delay line');
end

% 0.0 delay is critical for some algorithms, TODO: handle this differently
fdt( fdt == 0 ) = 1E-6;  

% There is no post processing stage if the Farrow Structure used
if ~strcmp( fracdelay.pre.method, 'farrow' )
  % === Post Processing ====================================================
  h = zeros(Nh, channels);  % fractional delay filter
  switch fracdelay.filter
    case 'zoh'
      % === Zero-Order-Hold (Integer Delays) ===============================
      idt = ceil(dt);  % round up to next integer delay
      fdt = 0.0;  % 
      h = ones(1, channels);
    case 'lagrange'
      % ==== Lagrange Polynomial Interpolator ==============================
      c = lagrange_polynoms(0:(Nh-1));  %
      for sdx=1:Nh
        h(sdx,:) = polyval(c(sdx,:),fdt);
      end
    case 'least_squares'
      % ==== Least Squares Interpolation Filter ============================
      for cdx=1:channels
        h(:,cdx) = general_least_squares(Nh,fdt(cdx),0.90);
      end
    otherwise
      disp('Delayline: Unknown Pre-Processing method for delay line');
  end
  
  for cdx=1:channels
    delayline(:,cdx) = filter(h(:,cdx),1,delayline(:,cdx));
  end  
end

%% ===== Integer Delay ================================================
% Handling of too long delay values (returns vector of zeros)
idt(abs(idt) > samples) = samples;

% Handle positive or negative delays
for cdx=1:channels
  if idt(cdx)>=0
    delayline(:,cdx) = [zeros(idt(cdx),1); weight(cdx)*delayline(1:end-idt(cdx),cdx)];
  else
    delayline(:,cdx) = [weight(cdx)*delayline(-idt(cdx)+1:end,cdx); zeros(-idt(cdx),1)];
  end
end
sig = delayline(1:rfactor:samples,:);

%%
% Undo reshaping [N M*C] => [M C N]
if reshaped
  sig = reshape(sig',[M C size(sig,1)]);
end