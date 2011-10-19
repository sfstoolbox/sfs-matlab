function [itd,idxleft,idxright] = extract_itd(insigleft,insigright,fs)
%EXTRACTITD Extract the ITD between the two given signals
%   Usage: itd = extractitd(insigleft,insigrighti,fs)
%
%   Input parameters:
%       insigleft   - left ear signal. This can also be a matrix containing
%                     left signals for different frequency bands
%       insigright  - the same as insigleft, but for the right ear
%       fs          - sampling rate
%
%   Output parameters:
%       itd         - ITD for the given signals. A single value for two
%                     given signals or a vector with values for every
%                     frequency band
%
%   EXTRACTITD(insigleft,insigright,fs) extractes the ITD between the left and
%   right signal(s) by using an edge detection algorithm to identify the
%   first non-zero entry in both IRs and then calculating the time
%   difference.
%
%R gaik1993, sandvad1994, lindau2010

% AUTHOR: Hagen Wierstorf


%% ------ Checking of input parameters -----------------------------------
nargmin = 3;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(insigleft,insigright);
isargpositivescalar(fs);
if size(insigright)~=size(insigright)
    error('%s: insigleft and insigright have to be the same size!', ...
        upper(mfilename));
end


%% ------ Computation ----------------------------------------------------

% Extract the envelope of the input signals
insigleft = abs(hilbert(insigleft));
insigright = abs(hilbert(insigright));

% See if we have more than one frequency channel in the insig
itd = zeros(1,size(insigleft,2));
idxleft = itd;
idxright = itd;
for ii = 1:size(insigleft,2)

    % Treshold after sandvad1994 (5% of maximum in each IR)
    % NOTE: I have changed it to 10%
    tresholdleft = 0.10 * max(insigleft(:,ii));
    tresholdright = 0.10 * max(insigright(:,ii));
    % Ten fold upsampling (after lindau2010) to have a smoother output
    resampleft = resample(insigleft(:,ii),10*fs,fs);
    resampright = resample(insigright(:,ii),10*fs,fs);

    % Find maximum and use the maximum to calculate the ITD
    %[maxleft idxleft(ii)] = max(insigleft(:,ii));
    %[maxright idxright(ii)] = max(insigright(:,ii));

    idxleft(ii) = find(resampleft > tresholdleft,1,'first');
    idxright(ii) = find(resampright > tresholdright,1,'first');

    % Calculate ITD
    itd(ii) = (idxleft(ii)-idxright(ii))/(10*fs);

end
