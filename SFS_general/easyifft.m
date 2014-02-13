function outsig = easyifft(amplitude,phase)
%EASYIFFT calculates the inverse FFT
%
%   Usage: outsig = easyifft(amplitude,phase)
%
%   Input parameters:
%       amplitude   - the amplitude spectrum
%       phase       - the phase spectrum / rad
%
%   Output parameters:
%       outsig      - a one channel signal
%
%   EASYIFFT(amplitude,phase) generates the corresponding waveform from the
%   amplitude and phase spectra using ifft.
%
%   see also: easyfft, ifft

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


%% ===== Checking input arguments ========================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
[amplitude,phase] = column_vector(amplitude,phase);


%% ===== Regenerating wave form from spectrum ============================
% Length of the signal to generate
samples = 2 * (length(amplitude)-1);

% Rescaling (see easyfft)
amplitude = amplitude/2 * samples;

% Mirror the amplitude spectrum
amplitude = [ amplitude; amplitude(end-1:-1:2) ];

% Mirror the phase spectrum and build the inverse (why?)
phase = [ phase; -1*phase(end-1:-1:2) ];

% Convert to complex spectrum
compspec = amplitude .* exp(1i*phase);

% Build the inverse fft and use only the real part
outsig = real( ifft(compspec) );
