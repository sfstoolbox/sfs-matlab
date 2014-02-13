function outsig = noise(samples,nsigs,noisetype)
%NOISE creates a white noise signal
%
%   Usage: outsig = noise(samples,[nsigs,[type]])
%
%   Input parameters:
%       samples     - length of the noise signal / samples
%       nsigs       - number of noise signals, default: 1
%       type        - type of noise:
%                         'white' (default)
%                         'pink'
%                         'red'
%                         'brown'
%
%   Output parameters:
%       outsig  - samples x nsigs white noise signal, default: nsigs = 1.
%
%   NOISE(samples,nsigs,type) generate a noise signal of type with a length of
%   samples and nsigs columns. The default value is nsigs = 1; the default value
%   for type is white.
%
%   NOTE: be sure you set the state of the random number generator in your
%   Matlab session, otherwise you will get the same noise every time you
%   restart Matlab.
%   The number generator is been set by:
%
%       randn('state',sum(100*clock));  % in Matlab < 7.7
%
%       s = RandStream.create('mt19937ar','seed',sum(100*clock));
%       RandStream.setDefaultStream(s);     % in Matlab >= 7.7
%
%   see also: click, irn

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
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


%% ===== Checking input parameters =======================================
nargmin = 1;
nargmax = 3;
narginchk(nargmin,nargmax);
if nargin<3
    noisetype = 'white';
elseif nargin<2
    nsigs = 1;
end
isargpositivescalar(samples,nsigs)
isargchar(noisetype)


%% ===== Computation =====================================================
%
switch noisetype
    case 'white'
        % Generate white noise
        outsig = randn(samples,nsigs);

        % Generate white noise using FFT
        % Amplitude = 1
        %amplitude = ones(samples,1);
        % Generate random phase
        %phase = 2*pi*rand(samples,1);
        % Generate noise signal using fft
        %outsig = easyifft(amplitude,phase);

    case 'pink'
        % Handle trivial condition
        if samples==1
            outsig = ones(1,nsigs);
            return;
        end
        fmax = floor(samples/2)-1;
        f = (2:(fmax+1)).';
        % 1/f amplitude factor
        a = 1./sqrt(f);
        % Random phase
        p = randn(fmax,nsigs) + 1i*randn(fmax,nsigs);
        sig = repmat(a,1,nsigs).*p;
        outsig = ifftreal([ones(1,nsigs); sig; ...
            1/(fmax+2)*ones(1,nsigs)],samples);

    case 'red'
        outsig = cumsum(randn(samples,nsigs));
        
    case 'brown'
        outsig = cumsum(randn(samples,nsigs));

    otherwise
        error('%s: unknown noise type %s.',upper(mfilename),noisetype)
end

% Scale output noise signal
outsig = norm_signal(outsig);
