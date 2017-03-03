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
%   See also: aliasing_frequency, secondary_source_positions

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


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax),
isargsecondarysource(x0);
if nargin==nargmax-1
    approx = false;
end


%% ===== Calculation ====================================================
% Calculate the distance to the nearest secondary source for all secondary
% sources
if size(x0,1)==1
    % If we have only one speaker
    dx0 = Inf;
else
    if approx
        % Approximate by using only the first <=100 sources
        x0 = x0(1:min(100,size(x0,1)),:);
    end
    dx0_single = zeros(size(x0,1),1);
    for ii=1:size(x0,1)
        % First secondary source position
        x01 = x0(ii,1:3);
        % All other positions
        x02 = [x0(1:ii-1,1:3); x0(ii+1:end,1:3)];
        % Get distance between x01 and all secondary sources within x02
        dist = vector_norm(bsxfun(@minus,x02,x01),2);
        % Get the smallest distance (which is the one to the next source)
        dx0_single(ii) = min(dist);
    end
    % Get the mean distance between all secondary sources
    dx0 = mean(dx0_single);
end
