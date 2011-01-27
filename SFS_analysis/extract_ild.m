function ild = extract_ild(insigleft,insigright)
%EXTRACTILD Extract the ILD between the two given signals
%   Usage: ild = extractild(insigleft,insigright)
%
%   Input parameters:
%       insigleft   - left ear signal. This can also be a matrix containing
%                     left signals for different frequency bands
%       insigright  - the same as insigleft, but for the right ear
%
%   Output parameters:
%       ild         - ILD for the given signals. A single value for two
%                     given signals or a vector with values for every
%                     frequency band
%
%   EXTRACTILD(insigleft,insigright) extractes the ILD between the left and
%   right signal(s) by subtracting the dB value of the left signal(s) from
%   the dB value of the right signal(s). 
%
%R gaik1993

% AUTHOR: Hagen Wierstorf


%% ------ Checking of input parameters -----------------------------------

error(nargchk(2,2,nargin));

if ~isnumeric(insigleft)
    error('%s: insigleft has to be a numeric signal!',upper(mfilename));
end
if ~isnumeric(insigright)
    error('%s: insigright has to be a numeric signal!',upper(mfilename));
end
if size(insigright)~=size(insigright)
    error('%s: insigleft and insigright have to be the same size!', ...
        upper(mfilename));
end


%% ------ Computation ----------------------------------------------------

% See if we have more than one frequency channel in the insig
ild = zeros(1,size(insigleft,2));
for ii = 1:size(insigleft,2)
    
    ild(ii) = rmsdb(insigright(:,ii))-rmsdb(insigleft(:,ii));
    
end
    
