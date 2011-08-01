function [amplitude,phase,f] = easyfft(insig,fs)
%EASYFFT Calculates the FFT of a signal
%   Usage: [amplitude,phase,f] = easyfft(insig,fs)
%
%   Input parameters:
%       insig       - one channel audio waveform
%       fs          - sampling rate (this value is used to calculate an 
%                     appropriate frequency axis)
%
%   Output parameters:
%       amplitude   - amplitude spectrum of the input signal
%       phase       - phase spectrum of the input signal
%       f           - corresponding frequency axis for the amplitude 
%                     spectrum (=> plot(f,amplitude)
%
%   EASYFFT(insig,fs) calculates the fft of the given signal.
%
%   see also: easyifft, fft
%

% AUTHOR: Hagen Wierstorf


%% ===== Check input arguments ===========================================
nargmin = 2;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargvector(insig)
if size(insig,2) > 1
   % Transpose
   insig = insig';
end;


%% ===== Calcualate spectrum =============================================
% Generate fast fourier transformation (=> complex output)
compspec = fft(insig);

% Length of the signal => number of points of fft
samples = length(insig);

% Get amplitude and phase spectra (and use only the first half of the
%>spectrum (Nyquist))
amplitude = abs(compspec(1:ceil(samples/2)));
phase = angle(compspec(1:ceil(samples/2)));

% Scale the amplitude (factor two, because we have cut off one half and 
%>divide by number of samples)
amplitude = 2*amplitude / samples;

% Calculate corresponding frequency axis
f = fs*(0:ceil(samples/2)-1)/samples;
