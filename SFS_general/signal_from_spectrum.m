function outsig = signal_from_spectrum(amplitude,phase,f,conf)
%SIGNAL_FROM_SPECTRUM creates time signal from single-sided spectrum
%
%   Usage: outsig = signal_from_spectrum(amplitude,phase,f,conf)
%
%   Input parameters:
%       amplitude   - the single-sided amplitude spectrum
%       phase       - the single-sided phase spectrum / rad
%       f           - the corresponding frequency vector
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       outsig      - one channel audio (time) signal
%
%   SIGNAL_FROM_SPECTRUM(amplitude,phase,f,conf) generates the time signal from
%   single-sided amplitude and phase spectra using ifft. It is the counterpart
%   of SPECTRUM_FROM_SIGNAL and not interchangeable with calling ifft.
%
%   See also: spectrum_from_signal, ifft

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


%% ===== Checking input arguments ========================================
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
[amplitude,phase] = column_vector(amplitude,phase);


%% ===== Configuration ===================================================
fs = conf.fs;


%% ===== Regenerating wave form from spectrum ============================
% Provided number of frequency bins
bins = length(f);

if f(end) == fs/2  % -> even time signal length
    % Length of the signal to generate
    samples = 2 * (bins-1);
    % Rescaling (see spectrum_from_signal())
    amplitude = [amplitude(1); amplitude(2:end-1)/2; amplitude(end)] * samples;
    % Mirror the amplitude spectrum ( 2*pi periodic [0, fs[ )
    amplitude = [amplitude; amplitude(end-1:-1:2)];
    % Mirror the phase spectrum and build the inverse (complex conjugate)
    phase = [phase; -1 * phase(end-1:-1:2)];

else  % -> odd time signal length
    % Length of the signal to generate
    samples = 2*bins - 1;
    % Rescaling (see signal_from_spectrum)
    amplitude = [amplitude(1); amplitude(2:end)/2] * samples;
    % Mirror the amplitude spectrum ( 2*pi periodic [0, fs-bin] )
    amplitude = [amplitude; amplitude(end:-1:2)];
    % Mirror the phase spectrum and build the inverse (complex conjugate)
    phase = [phase; -1*phase(end:-1:2)];
end

% Convert to complex spectrum
compspec = amplitude .* exp(1i*phase);

% Build the inverse fft and assume spectrum is conjugate symmetric
outsig = real(ifft(compspec));
