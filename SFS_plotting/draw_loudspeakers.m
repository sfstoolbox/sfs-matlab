function draw_loudspeakers(x0,win,conf)
%DRAW_LOUDSPEAKERS draws loudspeaker symbols or "x" at the given positions
%
%   Usage: draw_loudspeakers(x0,[win,[conf]])
%
%   Input options:
%       x0          - positions and directions of the loudspeakers / m
%       win         - tapering window, which is the activity of the loudspeaker
%                     (default: 1)
%       conf        - optional configuration struct (see SFS_config)
%
%   DRAW_LOUDSPEAKERS(x0,win) draws loudspeaker symbols or crosses at the given
%   secondary source positions. This can be controlled by the
%   conf.plot.realloudspeakers setting. The loudspeaker symbols are pointing in
%   their given direction. In addition to the secondary source positions, the
%   activity of the single secondary sources can be given by the win vector. For
%   every secondary source it can contain a value between 0 and 1. 1 is fully
%   active. If only one value if given, it is used for all secondary sources.
%
%   see also: plot_wavefield

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input parameter =====================================
nargmin = 1;
nargmax = 3;
narginchk(nargmin,nargmax);
isargsecondarysource(x0)
nls = size(x0,1);
if nargin<nargmax-1
    win = ones(nls,1);
elseif length(win)==1
    win = win*ones(nls,1);
end
isargvector(win);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
p.realloudspeakers = conf.plot.realloudspeakers;
p.lssize = conf.plot.lssize;


%% ===== Plotting ========================================================
% Plot only "x" at the loudspeaker positions
if(~p.realloudspeakers)
    plot(x0(:,1),x0(:,2),'wx','linewidth',1,...
        'Color',[.01 .01 .01]);
else

    w = p.lssize;
    h = p.lssize;

    % Plot a real speaker symbol at the desired position
    % vertex coordinates
    %v1 = [-h/2 -h/2 0 0 -h/2 ; -w/2 w/2 w/2 -w/2 -w/2];
    %v2 = [0 h/2 h/2 0 ; -w/6 -w/2 w/2 w/6];
    v1 = [-h -h -h/2 -h/2 -h ; -w/2 w/2 w/2 -w/2 -w/2];
    v2 = [-h/2 0 0 -h/2 ; -w/6 -w/2 w/2 w/6];

    hold on;

    % draw loudspeakers
    for n=1:nls

        % Get the azimuth direction of the secondary sources
        [phi,r_tmp] = cart2pol(x0(n,4),x0(n,5));

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

        if(win(n)>0)
            % Set fill color for active loudspeakers
            % Scale the color. sc = 1 => black. sc = 0.5 => gray.
            sc = 0.5;
            fc = [(1-sc*win(n)), ...
                  (1-sc*win(n)), ...
                  (1-sc*win(n))];
        else
            % set fill color to white for inactive loudspeakers
            fc = [1,1,1];
        end

        % Draw speakers
        fill(v01(1,:),v01(2,:),fc);
        fill(v02(1,:),v02(2,:),fc);

    end

    hold off;
    
    % set equal axis ratio
    axis image

end
