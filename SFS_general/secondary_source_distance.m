function [dx0,dx0_single] = secondary_source_distance(x0,approx)
%SECONDARY_SOURCE_DISTANCE calculates the average distance between the secondary
%sources
%
%   Usage: dx0 = secondary_source_distance(x0,[approx])
%
%   Input parameters:
%       x0          - secondary sources / m
%       approx      - flag (0/1). If set only the first 100 secondary sources
%                     are used to calculate the distance. This can speed up a
%                     lot the calculation if you are using many secondary
%                     sources. Default: 0.
%
%   Output parameters:
%       dx0         - average distance of secondary sources / m
%       dx0_single  - vector containing the minimum distances of
%                     all secondary sources to the other ones / m
%
%   SECONDARAY_SOURCE_DISTANCE,approx(x0) calculates the average distance dx0
%   between the given secondary sources. First, the distance to its nearest
%   source is calculated for every single secondary source, afterwads the mean
%   about these values is returned.
%
%   see also: aliasing_frequency, secondary_source_positions

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
narginchk(nargmin,nargmax),
isargsecondarysource(x0);
if nargin==nargmax-1
    approx = false;
end


%% ===== Calculation ====================================================
% calculate the distance to the nearest secondary source for all secondary
% sources
if size(x0,1)==1
    % if we have only one speaker
    dx0 = Inf;
else
    if approx
        % approximate by using only the first <=100 sources
        x0 = x0(1:min(100,size(x0,1)),:);
    end
    dx0_single = zeros(size(x0,1),1);
    for ii=1:size(x0,1)
        % first secondary source position
        x01 = x0(ii,1:3);
        % all other positions
        x02 = [x0(1:ii-1,1:3); x0(ii+1:end,1:3)];
        % get distance between x01 and all secondary sources within x02
        dist = vector_norm(bsxfun(@minus,x02,x01),2);
        % get the smallest distance (which is the one to the next source)
        dx0_single(ii) = min(dist);
    end
    % get the mean distance between all secondary sources
    dx0 = mean(dx0_single);
end
