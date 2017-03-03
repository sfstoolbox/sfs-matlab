function draw_loudspeakers(x0,dimensions,conf)
%DRAW_LOUDSPEAKERS draws loudspeaker symbols or "o" at the given positions
%
%   Usage: draw_loudspeakers(x0,[dimensions],conf)
%
%   Input options:
%       x0          - positions and directions of the loudspeakers / m
%       dimensions  - dimension defining the plane in which the loudspeaker
%                     symbol should be plotted. For example [1 1 0] corresponds
%                     to the xy-plane. Default: [1 1 0]
%       conf        - configuration struct (see SFS_config)
%
%   DRAW_LOUDSPEAKERS(x0,dimensions,conf) draws loudspeaker symbols or filled
%   points at the given secondary source positions. This can be controlled by
%   the conf.plot.realloudspeakers setting. The loudspeaker symbols are pointing
%   in their given direction and are colored accordingly to their weights,
%   whereby the weights are scaled to an absolute maximum of 1 before.
%
%   See also: plot_sound_field

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


%% ===== Checking of input parameter =====================================
nargmin = 2;
nargmax = 3;
narginchk(nargmin,nargmax);
isargsecondarysource(x0)
nls = size(x0,1);
if nargin<nargmax
    conf = dimensions;
    dimensions = [1 1 0];
end


%% ===== Configuration ===================================================
p.realloudspeakers = conf.plot.realloudspeakers;
p.lssize = conf.plot.lssize;


%% ===== Plotting ========================================================
% Check in which plane we should plot the secondary sources and shift them
% around accordingly
if ~any(dimensions)
    return;
elseif ~dimensions(1)
    x0(:,1:2) = x0(:,2:3);
elseif ~dimensions(2)
    x0(:,2) = x0(:,3);
end

% Plot only "o" at the loudspeaker positions, use this as default for all cases
% that are not the x-y-plane
if ~p.realloudspeakers || ~(dimensions(1)&&dimensions(2))

    if p.realloudspeakers && ~(dimensions(1)&&dimensions(2))
        warning('%s: Real loudspeaker can only be drawn in the x-y-plane', ...
            upper(mfilename));
    end
    % Plot all secondary sources with the same color + symbol
    plot(x0(:,1),x0(:,2),'o', ...
        'MarkerFaceColor','k', ...
        'MarkerEdgeColor','k', ...
        'MarkerSize',4);

else

    % Weightings of the single loudspeakers (scale maximum to 1)
    weights = x0(:,7) / max(abs(x0(:,7)));
    % Set fill color for active loudspeakers
    % Scale the color: sc = 1 => black, sc = 0.5 => gray.
    sc = 0.6;
    fc = [(1-sc.*weights), ...
          (1-sc.*weights), ...
          (1-sc.*weights)];

    % Plot a real speaker symbol at the desired position
    % vertex coordinates with height and width after lssize
    w = p.lssize;
    h = p.lssize;
    v1 = [-h -h -h/2 -h/2 -h ; -w/2 w/2 w/2 -w/2 -w/2 ; 0 0 0 0 0];
    v2 = [-h/2 0 0 -h/2 ; -w/6 -w/2 w/2 w/6 ; 0 0 0 0 ];

    hold on;

    % Draw loudspeakers
    for n=1:nls

        % Get the azimuth direction of the secondary sources
        [phi,~] = cart2pol(x0(n,4),x0(n,5));

        % Rotation matrix (orientation of the speakers)
        % R = [cos(phi(n)) -sin(phi(n));sin(phi(n)) cos(phi(n))];
        R = rotation_matrix(phi);
        for k=1:length(v1)
            vr1(:,k) = R * v1(:,k);
        end
        for k=1:length(v2)
            vr2(:,k) = R * v2(:,k);
        end

        % Shift
        v01(1,:) = vr1(1,:) + x0(n,1);
        v01(2,:) = vr1(2,:) + x0(n,2);
        v02(1,:) = vr2(1,:) + x0(n,1);
        v02(2,:) = vr2(2,:) + x0(n,2);

        % Draw loudspeakers
        fill(v01(1,:),v01(2,:),fc(n,:));
        fill(v02(1,:),v02(2,:),fc(n,:));

    end

    hold off;

end

% Set equal axis ratio
axis image
