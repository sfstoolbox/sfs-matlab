function [amplitude,phase,f] = easyfft(sig,conf)
%EASYFFT calculates the FFT of a signal and returns the corresponding frequency
%   axis
%
%   Usage: [amplitude,phase,f] = easyfft(sig,conf)
%          [amplitude,phase,f] = easyfft(sig)
%
%   Input parameters:
%       sig         - one channel audio waveform
%       conf        - optional struct containing configuration variables (see
%                     SFS_config for default values)
%
%   Output parameters:
%       amplitude   - amplitude spectrum of the input signal
%       phase       - phase spectrum of the input signal
%       f           - corresponding frequency axis for the amplitude 
%                     spectrum (=> plot(f,amplitude) (Hz)
%
%   EASYFFT(sig) calculates the amplitude and phase of the sig spectrum by using
%   the fast Fourier transformation. In addition to the amplitude and phase, the
%   corresponding frequency axis for a plot is returned.
%
%   see also: easyifft, fft
%

% AUTHOR: Hagen Wierstorf


%% ===== Check input arguments ===========================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
sig = column_vector(sig);

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


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
f = fs*(0:ceil(samples/2)-1)/samples;
