function varargout = easyfft(sig,conf)
%EASYFFT calculates the FFT of a signal and returns the corresponding frequency
%   axis
%
%   Usage: [amplitude,phase,f] = easyfft(sig,conf)
%
%   Input parameters:
%       sig         - one channel audio waveform
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       amplitude   - amplitude spectrum of the input signal
%       phase       - phase spectrum of the input signal / rad
%       f           - corresponding frequency axis for the amplitude
%                     spectrum (=> plot(f,amplitude) / Hz
%
%   EASYFFT(sig,conf) calculates the amplitude and phase of the sig spectrum by
%   using the fast Fourier transformation. In addition to the amplitude and
%   phase, the corresponding frequency axis for a plot is returned.
%
%   See also: easyifft, fft

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


%% ===== Check input arguments ===========================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
sig = column_vector(sig);
isargstruct(conf);


%% ===== Configuration ===================================================
fs = conf.fs;
useplot = conf.plot.useplot;


%% ===== Calcualate spectrum =============================================
% Generate fast fourier transformation (=> complex output)
compspec = fft(sig);

% Length of the signal => number of points of fft
samples = length(sig);

% Get amplitude and phase spectra (and use only the first half of the
%>spectrum (Nyquist))
amplitude = abs(compspec(1:ceil(samples/2)));
phase = angle(compspec(1:ceil(samples/2)));

% Scale the amplitude (factor two, because we have cut off one half and
%>divide by number of samples)
amplitude = 2*amplitude / samples;

% Calculate corresponding frequency axis
f = fs*(0:ceil(samples/2)-1)'/samples;

% Return values
if nargout>0, varargout{1}=amplitude; end
if nargout>1, varargout{2}=phase; end
if nargout>2, varargout{3}=f; end


%% ===== Plotting ========================================================
if nargout==0 || useplot
    figure; semilogx(f,20*log10(abs(amplitude)));
end
