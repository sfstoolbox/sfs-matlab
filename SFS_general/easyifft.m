function outsig = easyifft(amplitude,phase)
%EASYIFFT calculates the inverse FFT
%
%   Usage: outsig = easyifft(amplitude,phase)
%
%   Input parameters:
%       amplitude   - the amplitude spectrum
%       phase       - the phase spectrum
%
%   Output parameters:
%       outsig      - a one channel signal
%
%   EASYIFFT(amplitude,phase) generates the corresponding waveform from the
%   amplitude and phase spectra using ifft.
%
%   see also: easyfft, ifft
%

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking input arguments ========================================
nargmin = 2;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
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
