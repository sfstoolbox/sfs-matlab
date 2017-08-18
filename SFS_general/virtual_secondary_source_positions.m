function xv = virtual_secondary_source_positions(x0,xs,src,conf)
%VIRTUAL_SECONDARY_SOURCE_POSITIONS generates the positions and directions of a
%   virtual secondary source distribution
%
%   Usage: xv = virtual_secondary_source_positions(x0,xs,src,conf)
%
%   Input options:
%       x0          - positions, directions and weights of real secondary
%                     sources [nx7]
%       xs          - position and for focused sources also direction of the
%                     desired source model / m [1x3] or [1x6]
%       src         - source type of the target field
%                       'pw' - plane wave (xs is the direction of the
%                              plane wave in this case)
%                       'ps' - point source
%                       'fs' - focused source (not supported, yet)
%       conf        - configuration struct (see SFS_config)
%
%   Output options:
%       xv          - virtual secondary source positions, directions and
%                     weights / m
%
%   VIRTUAL_SECONDARY_SOURCE_POSITIONS(x0,xs,src,conf) generates the positions
%   and directions xv of virtual secondary sources for a local area geometry
%   (conf.localsfs.geometry) and local area size (conf.localsfs.size).
%   The direction of the virtual sources is given as their unit vectors
%   pointing in the given direction.
%   Optionally (conf.localwfs_vss.consider_target_field == true), the
%   algorithm takes the sound field, which is to be reproduced, into account
%   for the positioning.
%   Optionally (conf.localwfs_vss.consider_secondary_sources == true), the
%   algorithm takes the positions of the real loudspeakers into account.
%
%   See also: driving_function_mono_localwfs, driving_function_imp_localwfs

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


%% ===== Checking of input  parameters ===================================
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargxs(xs);
isargsecondarysource(x0);
isargchar(src);
isargstruct(conf);


%% ===== Configuration ===================================================
virtualconf = conf;
virtualconf.secondary_sources.size = conf.localwfs_vss.size;
virtualconf.secondary_sources.center = conf.localwfs_vss.center;
virtualconf.secondary_sources.geometry = conf.localwfs_vss.geometry;
virtualconf.secondary_sources.number = conf.localwfs_vss.number;
virtualconf.secondary_sources.grid = conf.localwfs_vss.grid;

geometry                    = conf.localwfs_vss.geometry;
nls                         = conf.localwfs_vss.number;
consider_secondary_sources  = conf.localwfs_vss.consider_secondary_sources;
consider_target_field       = conf.localwfs_vss.consider_target_field;


%% ===== Main ============================================================

if consider_target_field || consider_secondary_sources
    % =====================================================================
    % Adaptive positioning of virtual secondary sources
    % =====================================================================

    Rl = virtualconf.secondary_sources.size/2;  % radius of local area
    xl = virtualconf.secondary_sources.center;  % center of local area

    % Determine vector poiting towards source
    if strcmp('pw',src)
        % === Plane wave ===
        ns = bsxfun(@rdivide,-xs,vector_norm(xs,2));
    elseif strcmp('ps',src) || strcmp('ls',src)
        % === Point source ===
        ns = bsxfun(@rdivide,xs-xl,vector_norm(xs-xl,2));
    elseif strcmp('fs',src)
        % === Focused source ===
        ns = bsxfun(@rdivide,xl-xs,vector_norm(xs-xl,2));
    end
    phis = atan2(ns(2),ns(1));  % azimuth angle of ns

    if strcmp('circle',geometry) || strcmp('circular',geometry)
        % =====================================================================
        % Virtual circular Array
        % =====================================================================
        % Valid arc of virtual secondary sources

        % CONSTRAINT 1 ========================================================
        % Consider the target field which shall reproduced
        phid = pi;
        if consider_target_field
            if strcmp('pw',src)
                phid = pi/2;
            elseif strcmp('fs', src)
                phid = acos(vector_norm(xs-xl,2)./Rl);
            else
                % 1/2 opening angle of cone spanned by local area and virtual source
                phid = acos(Rl./vector_norm(xs-xl,2));
            end
        end

        % CONSTRAINT 2 ======================================================
        % Consider the position of the real loudspeaker array
        delta_max = phid;
        delta_min = -phid;
        if (consider_secondary_sources && ~isempty(x0))
            delta_max = 0;
            delta_min = 0;
            % For each secondary source
            for idx=1:size(x0,1)
                xc0 = x0(idx,1:3) - xl;  % vector from secondary source to local area
                Rc0 = vector_norm(xc0,2);  % distance from secondary source to local area
                nc0 = xc0./Rc0;  % normal vector from secondary source to local area
                % 1/2 opening angle of cone spanned by local area and secondary source
                phix0 = acos(Rl./Rc0);

                phiso = asin(ns(1)*nc0(2) - ns(2)*nc0(1));  % angle between ns and nc0
                delta_max = max(delta_max, phiso + phix0);
                delta_min = min(delta_min, phiso - phix0);
            end
            delta_max = min(delta_max, phid);
            delta_min = max(delta_min, -phid);
        end

        delta_offset = eps;

        % SOURCE POSITIONING ==================================================
        % === Equi-angular sampling on valid arc ===
        phi = phis + linspace(delta_min + delta_offset,delta_max-delta_offset, nls).';

        % Elevation angles
        theta = zeros(nls,1);
        % Positions of the secondary sources
        [cx,cy,cz] = sph2cart(phi,theta,Rl);
        xv(:,1:3) = [cx,cy,cz] + repmat(xl,nls,1);
        % Direction of the secondary sources
        xv(:,4:6) = direction_vector(xv(:,1:3),repmat(xl,nls,1).*ones(nls,3));
        % Equal weights for all sources
        xv(:,7) = ones(nls,1) * (delta_max-delta_min)*Rl/(nls-1);

    elseif strcmp('linear',geometry)
        % =====================================================================
        % Virtual linear Array
        % =====================================================================
        % valid line of virtual secondary sources

        % CONSTRAINT 1 ========================================================
        % Consider the target field which shall reproduced
        if consider_target_field
            nd = ns;
            ndorth = ns*[0 1 0; -1 0 0; 0 0 1];
        else
            nd = [0, 1, 0];
            ndorth = [-1, 0, 0];
        end

        % CONSTRAINT 2 ========================================================
        % Consider the position of the real loudspeaker array
        if (consider_secondary_sources && ~isempty(x0))
            xmax = inf;
            xmin = -inf;
            % For each secondary source
            for idx=1:size(x0,1)
                % This calculates the parameter t for the intersection between the linear
                % virtual array and the plane spanned by one loudspeaker (position +
                % orientation)
                n0 = x0(idx, 4:6);
                ntmp = (ndorth*n0');
                if ntmp ~= 0
                    t = (x0(idx,1:3) - xl)*n0'./ntmp;
                    if t > 0
                        xmax = min(xmax, t);
                    else
                        xmin = max(xmin, t);
                    end
                end
            end
            xmax = min(xmax, Rl);
            xmin = max(xmin, -Rl);
        else
            xmax = Rl;
            xmin = -Rl;
        end

        % SOURCE POSITIONING ==================================================
        % === Equi-distant sampling on valid line ===
        x = linspace(xmin, xmax, nls);

        % Positions of the secondary sources
        xv(:,1:3) = repmat(xl, nls, 1) + x'*ndorth;
        % Direction of the secondary sources
        xv(:,4:6) = repmat(-nd, nls, 1);
        % Equal weights for all sources
        xv(:,7) = ones(nls,1)*(xmax-xmin)/(nls-1);
    else
        xv = secondary_source_positions(virtualconf);
        if consider_target_field
            xv = secondary_source_selection(xv, xs, src);
        end
    end
else
    % =====================================================================
    % Non-adaptive positioning of  virtual secondary sources
    % =====================================================================
    % Just position the virtual secondary sources as you would do it with the
    % loudspeakers. This ignores the virtual source position, the size of
    % the local listening area and the position of the real loudspeakers
    xv = secondary_source_positions(virtualconf);
end
