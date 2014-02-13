function win = tapering_window(x0,conf)
%TAPWIN generate a tapering window for a loudspeaker array
%
%   Usage: win = tapering_window(x0,[conf])
%
%   Input parameters:
%       x0          - secondary sources / m
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       win     - tapering window (nlsx1)
%
%   TAPERING_WINDOW(xo,conf) generates a tapering window for a secondary source
%   distribution given by x0. The window is created from a squared Hann window.
%   If the secondary source distribution has some gaps, every joint part gets
%   its own tapering.
%   The mean distance of the secondary sources is calculated within this
%   function in order to identify edges of the array.
%
%   see also: secondary_source_position, sound_field_mono_wfs, hann

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


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
isargsecondarysource(x0);
if nargin==nargmax-1
    conf = SFS_config;
end
isargstruct(conf);


%% ===== Configuration ==================================================
usetapwin = conf.usetapwin;
tapwinlen = conf.tapwinlen;
geometry = conf.secondary_sources.geometry;


%% ===== Calculation =====================================================
% number of speakers
nls = size(x0,1);
% FIXME: at the moment the tapering window is not working for spherical arrays,
% because we are not able to find the edges of the array.
if usetapwin && nls>2 && ...
   ~(strcmp('sphere',geometry)||strcmp('spherical',geometry))
    win = ones(1,nls);
    % get the mean distance between secondary sources and the smallest distance
    % to neighbour source for every secondary source. Due to long computing time
    % for really large secondary source numbers, the distance is approximated by
    % the first 100 (or less) sources. If you don't want this behavior, change
    % the following command to dx0 = secondary_source_distance(x0,0);
    dx0 = secondary_source_distance(x0,1);
    % use only positions
    x0 = x0(:,1:3);
    % find the edges of the array
    edges = [];
    for ii=1:nls-1
        if norm(x0(ii,:)-x0(ii+1,:))>2*dx0
            edges = [edges; ii; ii+1];
        end
    end
    if norm(x0(end,:)-x0(1,:))>2*dx0
        edges = [edges; nls; 1];
    end
    % if we have any edges in our array apply a tapering window for every array
    % part, consiting of two edges
    if ~isempty(edges)
        % generate tapwin for every array part within the x0 vector
        for ii=2:length(edges)-2
            part_nls = edges(ii+1)-edges(ii)+1;
            win(edges(ii):edges(ii+1)) = part_hann_win(part_nls,tapwinlen);
        end
        % generate tapwin for every array part consiting of the first and the
        % last edge within the x0 vector
        if edges(1)==nls
            win = part_hann_win(nls,tapwinlen);
        else
            part_nls = edges(1) + nls-edges(end)+1;
            part_win = part_hann_win(part_nls,tapwinlen);
            win(1:edges(1)) = part_win(end-edges(1)+1:end);
            win(edges(end):end) = part_win(1:end-edges(end)+1);
        end
    end
else
    % If you want to use no tapering window:
    win = ones(1,nls);
end

% Ensure the window will be a column vector
win = column_vector(win);

end

function [win] = part_hann_win(nls,tapwinlen)
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
        win = hann_window(floor((lenwin-2)/2),floor((lenwin-2)/2),nls).^2;

    end
end
