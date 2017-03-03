function varargout = spectrum_from_signal(signal,conf)
%SPECTRUM_FROM_SIGNAL returns single-sided amplitude and phase spectra of signal
%
%   Usage: [amplitude,phase,f] = spectrum_from_signal(sig,conf)
%
%   Input parameters:
%       signal      - one channel audio (time) signal
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       amplitude   - single-sided amplitude spectrum of the input signal
%       phase       - single-sided phase spectrum of the input signal / rad
%       f           - corresponding frequency axis for the spectrum
%                     => plot(f,amplitude) / Hz
%
%   SPECTRUM_FROM_SIGNAL(signal,conf) calculates the single-sided amplitude and
%   phase spectra of a time signal by using the fast Fourier transformation.
%   In addition to the amplitude and phase, the corresponding frequency axis is
%   returned.
%
%   See also: signal_from_spectrum, fft

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


%% ===== Check input arguments ===========================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
signal = column_vector(signal);
isargstruct(conf);


%% ===== Configuration ===================================================
fs = conf.fs;
useplot = conf.plot.useplot;


%% ===== Calcualate spectrum =============================================
% Generate fast fourier transformation (=> complex output)
compspec = fft(signal);

% Length of the signal => number of points of fft
bins = length(signal);

if mod(bins,2)  % For odd signal length
    % Calculate corresponding frequency axis
    f = fs/bins * (0:(bins-1)/2)';
    % Get amplitude and phase spectra and use only the first half of the
    %>spectrum [0, fs/2[
    amplitude = abs(compspec(1:length(f)));
    phase = angle(compspec(1:length(f)));
    % Scale the amplitude (factor two for mirrored frequencies
    %>divide by number of bins)
    amplitude = [amplitude(1); 2*amplitude(2:end)] / bins;

else  % For even signal length
    % Calculate corresponding frequency axis
    f = fs/bins * (0:bins / 2)';
    % Get amplitude and phase spectra and use only the first half of the
    %>spectrum [0, fs/2]
    amplitude = abs(compspec(1:length(f)));
    phase = angle(compspec(1:length(f)));
    % Scale the amplitude (factor two for mirrored frequencies
    %>divide by number of bins)
    amplitude = [amplitude(1); 2*amplitude(2:end-1); amplitude(end)] / bins;
end

% Return values
if nargout>0, varargout{1}=amplitude; end
if nargout>1, varargout{2}=phase; end
if nargout>2, varargout{3}=f; end


%% ===== Plotting ========================================================
if nargout==0 || useplot
    figure; title('Spectrum');
    subplot(2,1,1)
    semilogx(f,20 * log10(abs(amplitude))); xlim([1, fs/2]);
    grid on; xlabel('frequency / Hz'); ylabel('amplitude / dB')
    subplot(2,1,2)
    semilogx(f,unwrap(phase)); xlim([1, fs/2]);
    grid on; xlabel('frequency / Hz'); ylabel('phase / rad')
end
