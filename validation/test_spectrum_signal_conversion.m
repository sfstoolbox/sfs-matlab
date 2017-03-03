function status = test_spectrum_signal_conversion(modus)
%TEST_SPECTRUM_SIGNAL_CONVERSION tests the spectrum to and from signal
%conversions
%
%   Usage: status = test_spectrum_signal_conversion(modus)
%
%   Input parameters:
%       modus    - 0: numerical
%                  1: visual
%
%   Output parameters:
%       status - true or false
%
%   test_spectrum_signal_conversion(modus) checks if spectrum_from_signal() and
%   signal_from_spectrum() works correctly.

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


status = false;


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);


%% ===== Configuration ===================================================
conf = SFS_config;
fs = conf.fs;
%% Create defined signal
t = 0 : 1/fs : 1 - 1/fs;  % 1 s
sin1 = sin(2*pi*50 * t);
sin2 = sin(2*pi*300 * t);
sin3 = sin(2*pi*1000 * t);
sin_sig = (sin1 + sin2 + sin3)';

even_sig = ones(8, 1);
odd_sig = ones(7, 1);
alias_sig = repmat([1; -1], fs/2, 1);

%% FFT
[sin_ampl, sin_phase, sin_f] = spectrum_from_signal(sin_sig, conf);
[even_ampl, even_phase, even_f] = spectrum_from_signal(even_sig, conf);
[odd_ampl, odd_phase, odd_f] = spectrum_from_signal(odd_sig, conf);
[alias_ampl, alias_phase, alias_f] = spectrum_from_signal(alias_sig, conf);

%% Check frequency bins
if modus
figure; semilogx(sin_f,20*log10(sin_ampl)); title('Sinus Mix')
figure; scatter(even_f, even_ampl); title('Even signal (f=0) FFT');
figure; scatter(odd_f, odd_ampl);title('Odd signal (f=0) FFT');
figure; scatter(alias_f, alias_ampl);title('Signal with f=fs/2 FFT');
end
%% IFFT
sin_outsig = signal_from_spectrum(sin_ampl, sin_phase, sin_f,conf);
even_outsig = signal_from_spectrum(even_ampl, even_phase, even_f,conf);
odd_outsig = signal_from_spectrum(odd_ampl, odd_phase, odd_f,conf);
alias_outsig = signal_from_spectrum(alias_ampl, alias_phase, alias_f,conf);

%% Check Output
sin_diff = sum(abs(sin_sig - sin_outsig));
even_diff = sum(abs(even_sig - even_outsig));
odd_diff = sum(abs(odd_sig - odd_outsig));
alias_diff = sum(abs(alias_sig - alias_outsig));

if sum(sin_diff + even_diff + odd_diff + alias_diff) < 10^(-8)
    status = true;
end

end

