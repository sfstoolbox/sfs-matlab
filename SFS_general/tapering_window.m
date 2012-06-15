function win = tapering_window(L,ls_activity,conf)
%TAPWIN generate a tapering window for a loudspeaker array
%
%   Usage: win = tapering_window(L,[ls_activity],[conf])
%
%   Input parameters:
%       L           - length of the loudspeaker array (m)
%       ls_activity - vector containing the activity of the loudspeakers from
%                     0..1
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       win     - tapering window (nlsx1)
%
%   TAPERING_WINDOW(L,ls_activity,conf) generates a tapering window for a linear
%   WFS loudspeaker array with a length of L. The window is created from a
%   squared Hann window. For circular arrays it is necessary to apply the
%   ls_activity option.
%
%   see also: wave_field_wfs_25d, ir_wfs_25d, secondary_source_selection, hann

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

% AUTHOR: Hagen Wierstorf, Sascha Spors
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$

% TODO: the array length is only used, if ls_activity is not given. Maybe there
% is a better way to handle this


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));
isargpositivescalar(L);
if nargin==nargmax-1
    % Check if the second argument is conf or ls_activity
    if isstruct(ls_activity)
        conf = ls_activity;
        ls_activity = ones(1,secondary_source_number(L,conf));
    else
        conf = SFS_config;
    end
elseif nargin==nargmax-2
    conf = SFS_config;
    if strcmp(conf.array,'circle')
        error(['%s: For circular arrays tapering_window needs the ls_activity. ', ...
            'If you have really all loudspeakers active use ', ...
            'conf.usetapwin=0.'],upper(mfilename));
    end
    % If no explicit loudspeaker activity is given mark all speakers as active
    ls_activity = ones(1,secondary_source_number(L,conf));
end
isargstruct(conf);
isargvector(ls_activity);


%% ===== Configuration ==================================================
usetapwin = conf.usetapwin;
tapwinlen = conf.tapwinlen;


%% ===== Calculation ====================================================
% Generate a squared Hann window and split it to the two edges of the array as
% shown below.
%    ------------------------------------------------------------
%   |                                                            |
% _|                                                              |_
%
% This procedure becomes more complicated if not all speakers are active and
% if we have a circular array or any other closed loudspeaker array. The code
% below can handle all cases, where the active loudspeaker have no gaps.
%
% Find active loudspeaker and create only a window for these loudspeakers
idx = (( ls_activity>0 ));
nls = length(ls_activity(idx));

if(usetapwin)
    % Length of window (given by the value of tapwinlen). The window will be
    % splitted to both sides of the loudspeaker array.
    lenwin = round(tapwinlen*nls)+2;
    %
    % Check if we have a to short window to apply it in a useful way. This can
    % be the case for very short loudspeaker arrays (as used in Wierstorf2010).
    if lenwin<4
        win = ones(1,nls);
    else
        % Create a squared Hann window with length lenwin
        %    -
        %   | |
        % _|   |_
        hannwin = hann(lenwin).^2;
        % Create tapering window
        % NOTE: the first and the last entry generated by hann are 0. Therefore we
        % are using the squared Hann window from its second sample until its
        % second last one.
        win = [hannwin(2:ceil(end/2))' ...
               ones(1,nls-lenwin+2) ...
               hannwin(ceil(end/2)+1:end-1)'];
    end
else
    % If you want to use no tapering window:
    win = ones(1,nls);
end

% If we have non active loudspeaker we have to move the tapering window to the
% right position. Also we have to check for closed arrays.
if length(ls_activity)~=length(ls_activity(idx))
    % Look for the first inactive and for the first active loudspeaker
    idx1 = find(ls_activity==0,1,'first');
    idx2 = find(ls_activity==1,1,'first');
    if idx1~=1
        % If the first loudspeaker is active we apply the window from here on
        % until the first inactive. If there were additional active loudspeakers
        % at the end we know that we have a close array and the rest of the
        % window is applied at the end in a way that the tapering window will be
        % correct.
        win = [win(end-idx1+2:end), ...
               zeros(1,length(ls_activity)-length(ls_activity(idx))), ...
               win(1:end-idx1+1)];
    else
        % If we have an inactive loudspeaker at the beginning place the window
        % in the middle.
        win = [zeros(1,length(1:idx2-1)), ...
               win, ...
               zeros(1,length(ls_activity)-length(win)-length(1:idx2-1))];
    end
end

% Ensure the window will be a column vector
win = column_vector(win);
