function draw_loudspeakers(x0,dimensions,conf)
%DRAW_LOUDSPEAKERS draws loudspeaker symbols or "x" at the given positions
%
%   Usage: draw_loudspeakers(x0,[dimensions],[conf])
%
%   Input options:
%       x0          - positions and directions of the loudspeakers / m
%       dimensions  - dimension defining the plane in which the loudspeaker
%                     symbol should be plotted. For example [1 1 0] corresponds
%                     to the xy-plane
%       conf        - optional configuration struct (see SFS_config)
%
%   DRAW_LOUDSPEAKERS(x0,dimensions) draws loudspeaker symbols or crosses at
%   the given secondary source positions. This can be controlled by the
%   conf.plot.realloudspeakers setting. The loudspeaker symbols are pointing in
%   their given direction.
%
%   see also: plot_wavefield

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


%% ===== Checking of input parameter =====================================
nargmin = 1;
nargmax = 3;
narginchk(nargmin,nargmax);
isargsecondarysource(x0)
nls = size(x0,1);
if nargin==nargmax-1
    if isstruct(dimensions)
        conf = dimensions;
        dimensions = [1 1 0];
    else
        conf = SFS_config;
    end
elseif nargin==nargmax-2
    conf = SFS_config;
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

% Weightings of the single sources
win = x0(:,7);

% Plot only "x" at the loudspeaker positions, use this as default for all cases
% that are not the x-y-plane
if ~p.realloudspeakers || ~(dimensions(1)&&dimensions(2))
    if p.realloudspeakers && ~(dimensions(1)&&dimensions(2))
        warning('%s: Real loudspeaker can only be drawn in the x-y-plane', ...
            upper(mfilename));
    end
    % plot all secondary sources with the same color + symbol
    plot(x0(:,1),x0(:,2),'o', ...
        'MarkerFaceColor','k', ...
        'MarkerEdgeColor','k', ...
        'MarkerSize',4);
else
    % Set fill color for active loudspeakers
    % Scale the color. sc = 1 => black. sc = 0.5 => gray.
    sc = 0.6;
    fc = [(1-sc.*win), ...
          (1-sc.*win), ...
          (1-sc.*win)];

    w = p.lssize;
    h = p.lssize;

    % Plot a real speaker symbol at the desired position
    % vertex coordinates
    %v1 = [-h/2 -h/2 0 0 -h/2 ; -w/2 w/2 w/2 -w/2 -w/2];
    %v2 = [0 h/2 h/2 0 ; -w/6 -w/2 w/2 w/6];
    v1 = [-h -h -h/2 -h/2 -h ; -w/2 w/2 w/2 -w/2 -w/2 ; 0 0 0 0 0];
    v2 = [-h/2 0 0 -h/2 ; -w/6 -w/2 w/2 w/6 ; 0 0 0 0 ];
    
    hold on;

    % draw loudspeakers
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

        % shift
        v01(1,:) = vr1(1,:) + x0(n,1);
        v01(2,:) = vr1(2,:) + x0(n,2);
        v02(1,:) = vr2(1,:) + x0(n,1);
        v02(2,:) = vr2(2,:) + x0(n,2);


        % Draw speakers
        fill(v01(1,:),v01(2,:),fc(n,:));
        fill(v02(1,:),v02(2,:),fc(n,:));

    end

    hold off;
    
end

% set equal axis ratio
axis image
