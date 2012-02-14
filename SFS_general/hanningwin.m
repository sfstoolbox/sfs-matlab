function win = hanningwin(onset,offset,nsamples)
%HANNINGWIN generates a hanning window with on and off ramp
%
%   Usage: win = hanningwin(onset,offset,nsamples)
%
%   Input parameters:
%       onset		- onset in samples (0 for no onset)
%       offset		- offset in samples (0 for no offset)
%       nsamples    - length of the whole window (including on- and offset)
%
%   Output parameters:
%       win         - a hanning window (nsamples x 1) for multiplication 
%                     with the desired signal to be windowed
%
%   see also: click
%

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking input parameters =======================================
nargmin = 3;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));
isargpositivescalar(onset,offset,nsamples)
if nsamples<onset+offset
    error('%s: nsamples has to be greater than onset+offset.',...
        upper(mfilename));
end


%% ===== Computation =====================================================
onset=ceil(onset);
offset=ceil(offset);
%
% The hanning window looks like this
%
%           ###
%         ##   ##
%        #       #
%       #         #
%       #         #
%      #           #
%     #             # 
%   ##               ##
% ##                   ##
%
% So we need only one half!
% Therefore we will only use the first half of onsetwin and the second
% half of offsetwin.
%
%% Generate onset window
if onset==0
    onsetwin = [];
else
    tmp = hanning(2*onset);
    onsetwin = tmp(1:onset);
end
% Generate offset window
if offset==0
    offsetwin = [];
else
    tmp = hanning(2*offset);
    offsetwin = tmp(offset+1:end);
end

% Generate the complete window
win = [ onsetwin; ...
        ones(nsamples-onset-offset,1); ...
        offsetwin ];
