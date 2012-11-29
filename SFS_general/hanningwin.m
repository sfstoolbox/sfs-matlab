function win = hanningwin(onset,offset,nsamples)
%HANNINGWIN generates a hanning window with on and off ramp
%
%   Usage: win = hanningwin(onset,offset,nsamples)
%
%   Input parameters:
%       onset       - onset in samples (0 for no onset)
%       offset      - offset in samples (0 for no offset)
%       nsamples    - length of the whole window (including on- and offset)
%
%   Output parameters:
%       win         - a hanning window (nsamples x 1) for multiplication
%                     with the desired signal to be windowed
%
%   see also: click

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
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
