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
%   see also: click, tapering_window

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
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
    % generate an uneven window, see issue #18
    tmp = hann(2*onset+1);
    % disregard the first entry, because its zero
    onsetwin = tmp(2:onset+1);
end
% Generate offset window
if offset==0
    offsetwin = [];
else
    tmp = hann(2*offset+1);
    % disregard the last entry, becaus its zero
    offsetwin = tmp(offset+1:end-1);
end

% Generate the complete window
win = [ onsetwin; ...
        ones(nsamples-onset-offset,1); ...
        offsetwin ];
