function win = hann_window(onset,offset,nsamples)
%HANN_WINDOW generates a Hann window with on and off ramp
%
%   Usage: win = hann_window(onset,offset,nsamples)
%
%   Input parameters:
%       onset       - onset / samples (0 for no onset)
%       offset      - offset / samples (0 for no offset)
%       nsamples    - length of the whole window (including on- and offset)
%
%   Output parameters:
%       win         - a Hann window (nsamples x 1) for multiplication
%                     with the desired signal to be windowed
%
%   See also: click, tapering_window

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


%% ===== Checking input parameters =======================================
nargmin = 3;
nargmax = 3;
narginchk(nargmin,nargmax);
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
    % Generate an uneven window, see issue #18
    tmp = hann(2*onset+1);
    % Disregard the first entry, because its zero
    onsetwin = tmp(2:onset+1);
end
% Generate offset window
if offset==0
    offsetwin = [];
else
    tmp = hann(2*offset+1);
    % Disregard the last entry, becaus its zero
    offsetwin = tmp(offset+1:end-1);
end

% Generate the complete window
win = [ onsetwin; ...
        ones(nsamples-onset-offset,1); ...
        offsetwin ];
