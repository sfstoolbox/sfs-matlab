function [] = draw_loudspeakers(x0,y0,phi,ls_activity,conf)
%DRAW_LOUDSPEAKERS draws loudspeaker symbols or "x" at the given positions
%
%   Usage: draw_loudspeakers(x0,y0,phi,ls_activity,conf)
%          draw_loudspeakers(x0,y0,phi,ls_activity)
%          draw_loudspeakers(x0,y0,phi)
%
%   Input options:
%       x0,y0       - positions of the loudspeakers (m)
%       phi         - directions of the loudspeaker (rad)
%       ls_activity - activity of the loudspeaker
%       conf        - optional struct containing configuration variables (see
%                     SFS_config for default values)
%
%   DRAW_LOUDSPEAKERS(x0,y0,phi,ls_activity) draws loudspeaker symbols at
%   the given secondary source positions. The loudspeaker symbols are pointing
%   in the direction given by phi.
%
% AUTHOR: Sascha Spors, Hagen Wierstorf

%% ===== Checking of input parameter =====================================
nargmin = 3;
nargmax = 5;
error(nargchk(nargmin,nargmax,nargin));
isargvector(x0,y0,phi)
nLS = length(phi);
if(nargin<nargmax-1)
    ls_activity = zeros(1,nLS);
elseif(length(ls_activity)==1)
    ls_activity = ls_activity*ones(1,nLS);
end
isargvector(ls_activity);
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
    plot(x0(ls_activity>0),y0(ls_activity>0),'wx','linewidth',1,...
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
    for n=1:nLS

        % Rotation matrix (orientation of the speakers)
        % R = [cos(phi(n)) -sin(phi(n));sin(phi(n)) cos(phi(n))];
        R = rotation_matrix(phi(n)+pi/2);

        for k=1:length(v1)
            vr1(:,k) = R * v1(:,k);
        end

        for k=1:length(v2)
            vr2(:,k) = R * v2(:,k);
        end

        % shift
        v01(1,:) = vr1(1,:) + x0(n);
        v01(2,:) = vr1(2,:) + y0(n);

        v02(1,:) = vr2(1,:) + x0(n);
        v02(2,:) = vr2(2,:) + y0(n);

        if(ls_activity(n)>0)
            % Scale the color. sc = 1 => black. sc = 0.5 0> gray.
            sc = 0.5;
            fc = [(1-sc*ls_activity(n)), ...
                  (1-sc*ls_activity(n)), ...
                  (1-sc*ls_activity(n))];
            %fc = [ls_activity(n) ls_activity(n) ls_activity(n)];

            fill(v01(1,:),v01(2,:),fc);
            fill(v02(1,:),v02(2,:),fc);
        else
            h=line(v01(1,:),v01(2,:));
            set(h,'Color','k');
            %set(h,'Color',[.99 .99 .99]);
            h=line(v02(1,:),v02(2,:));
            set(h,'Color','k');
            %set(h,'Color',[.99 .99 .99]);
        end

    end

    hold off;

end
