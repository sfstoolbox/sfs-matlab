function xv = virtual_secondary_source_positions(x0,xs,src,conf)
%VIRTUAL_SECONDARY_SOURCE_POSITIONS Generates the positions and directions of a
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
%       conf        - optional configuration struct (see SFS_config)
%
%   Output options:
%       xv          - virtual secondary source positions, directions and
%                     weights / m
%
%   VIRTUAL_SECONDARY_SOURCE_POSITIONS(x0,xs,src,conf) generates the positions
%   and directions xv of virtual secondary sources for a local area geometry
%   (conf.localsfs.geometry) and local area size(conf.localsfs.size).
%   The direction of the virtual sources is given as their unit vectors
%   pointing in the given direction.
%   Optionally (conf.localsfs.vss.consider_target_field == true), the algorithm
%   takes the sound field, which is to be reproduced, into account for the
%   positioning.
%   Optionally (conf.localsfs.vss.consider_secondary_sources == true), the
%   algorithm takes the positions of the real loudspeakers into account. 

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

%% ===== Checking of input  parameters ===================================
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargxs(xs);
if ~isempty(x0)
  isargsecondarysource(x0);
end
isargchar(src);
if nargin<nargmax
  conf = SFS_config;
else
  isargstruct(conf);
end

%% ===== Configuration ===================================================
virtualconf = conf;
virtualconf.secondary_sources.size     = conf.localsfs.size;
virtualconf.secondary_sources.center   = conf.localsfs.center;
virtualconf.secondary_sources.geometry = conf.localsfs.vss.geometry;
virtualconf.secondary_sources.number   = conf.localsfs.vss.number;

geometry                    = conf.localsfs.vss.geometry;
sampling                    = conf.localsfs.vss.sampling;
logratio                    = conf.localsfs.vss.logratio;
nls                         = conf.localsfs.vss.number;
consider_local_area         = conf.localsfs.vss.consider_local_area;
consider_secondary_sources  = conf.localsfs.vss.consider_secondary_sources;
consider_target_field       = conf.localsfs.vss.consider_target_field;

%% ===== Main ============================================================

if consider_target_field || consider_secondary_sources || consider_local_area
  % =====================================================================
  % adaptive positioning of virtual secondary sources
  % =====================================================================

  Rl = virtualconf.secondary_sources.size/2;  % radius of local area
  xl = virtualconf.secondary_sources.center;  % center of local area

  % determine vector poiting towards source
  if strcmp('pw',src)
    % === Plane wave ===
    ns = bsxfun(@rdivide,-xs,vector_norm(xs,2));
  elseif strcmp('ps',src) || strcmp('ls',src)
    % === Point source ===
    ns = bsxfun(@rdivide,xs-xl,vector_norm(xs-xl,2));
  elseif strcmp('fs',src)
    % === Focused source ===
    to_be_implemented('focussed sources for virtual_secondary_source_positions');
  end
  phis = atan2(ns(2),ns(1));  % azimuth angle of ns

  if strcmp('circle',geometry) || strcmp('circular',geometry)
    % =====================================================================
    % virtual circular Array
    % =====================================================================
    % valid arc of virtual secondary sources

    % CONSTRAINT 1 ========================================================
    % consider the target field which shall reproduced
    phid = pi;
    if consider_target_field
      if strcmp('pw',src)
        phid = pi/2;
      else
        % 1/2 opening angle of cone spanned by local area and virtual source
        phid = acos(Rl./vector_norm(xs-xl,2));
      end
    end

    % CONSTRAINT 2 ======================================================
    % consider the position of the real loudspeaker array
    delta_max = phid;
    delta_min = -phid;
    if (consider_secondary_sources && ~isempty(x0))
      delta_max = 0;
      delta_min = 0;
      % for each secondary source
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
    switch (sampling)
      case 'equi'
        % === equi-angular sampling on valid arc ===
        phi = phis + linspace(delta_min + delta_offset,delta_max-delta_offset, nls).';
      case 'log'
        phi = log_spacing(delta_min + delta_offset, delta_max - delta_offset, 0, nls, logratio);
        phi = phis + phi;
      case 'scalar'
        x = linspace(sin(delta_min + delta_offset),sin(delta_max - delta_offset),nls).';
        phi = phis + asin(x);
      %case 'projective'
      %x = linspace(-3,3,nls).';
      %phi = phis + atan(x);
      otherwise
        error('%s: %s is not a supported sampling method for circular position!',upper(mfilename),method);
    end

    % Elevation angles
    theta = zeros(nls,1);
    % Positions of the secondary sources
    [cx,cy,cz] = sph2cart(phi,theta,Rl);
    xv(:,1:3) = [cx,cy,cz] + repmat(xl,nls,1);
    % Direction of the secondary sources
    xv(:,4:6) = direction_vector(xv(:,1:3),repmat(xl,nls,1).*ones(nls,3));
    % equal weights for all sources
    xv(:,7) = ones(nls,1);

  elseif strcmp('linear',geometry)
    % =====================================================================
    % virtual linear Array
    % =====================================================================
    % valid line of virtual secondary sources

    % CONSTRAINT 1 ========================================================
    % consider the target field which shall reproduced
    if consider_target_field
      nd = ns;
      ndorth = ns*[0 1 0; -1 0 0; 0 0 1];
      if consider_local_area
        if strcmp('pw',src)
          xd = Rl;
        else
          % 1/2 opening angle of cone spanned by local area and virtual source
          phid = acos(Rl./vector_norm(xs-xl,2));
          xd = (vector_norm(xs-xl,2)-Rl)./tan(phid);
        end
      else
        xd = Rl;
      end
    else
      xd = Rl;
      nd = [0, 1, 0];
      ndorth = [-1, 0, 0]; 
    end

    if consider_local_area
      Rd = Rl;
    else
      Rd = 0;
    end
    
    % CONSTRAINT 2 ========================================================
    % consider the position of the real loudspeaker array
    if (consider_secondary_sources && ~isempty(x0))
      xmax = inf;
      xmin = -inf;
      % for each secondary source
      for idx=1:size(x0,1)
        % this calculates the parameter t for the intersection between the linear
        % virtual array and the plane spanned by one loudspeaker (position +
        % orientation)
        n0 = x0(idx, 4:6);
        ntmp = (ndorth*n0');
        if ntmp ~= 0
          t = (x0(idx,1:3) - (xl + Rd*nd))*n0'./ntmp;
          if t > 0
            xmax = min(xmax, t);
          else
            xmin = max(xmin, t);
          end
        end
      end
      xmax = min(xmax, xd);
      xmin = max(xmin, -xd);
    else
      xmax = xd;
      xmin = -xd;
    end

    % SOURCE POSITIONING ==================================================
    switch (sampling)
      case 'equi'
        % === equi-distant sampling on valid line ===
        x = linspace(xmin, xmax, nls);
      case 'log'
        x = log_spacing(xmin, xmax, 0, nls, logratio);
      otherwise
        error('%s: %s is not a supported sampling method for linear sampling!',upper(mfilename),method);
    end

    % Positions of the secondary sources
    xv(:,1:3) = repmat(xl, nls, 1) + repmat(Rd*nd, nls, 1) + x'*ndorth;
    % Direction of the secondary sources
    xv(:,4:6) = repmat(-nd, nls, 1);
    % equal weights for all sources
    xv(:,7) = ones(nls,1);
  else
    error('%s: %s is not a supported positioning method!',upper(mfilename),method);
  end
else
  % =====================================================================
  % non-adaptive positioning of  virtual secondary sources
  % =====================================================================
  % just position the virtual secondary sources as you would do it with the
  % loudspeakers. This ignores the virtual source position, the size of
  % the local listening area and the position of the real loudspeakers

  xv = secondary_source_positions(virtualconf);
end

end

%% ===== Additional Helper Functions ====================================
function x = log_spacing(xmin, xmax, xcenter, N, ratio)

% variable indicating if number of samples is even/odd
even = mod(N,2) == 0;

if even
  N = N./2;
  s0 = -0.5;
else
  N = (N-1)./2;
  s0 = 0.0;
end

% variable indicating
xmax = abs(xmax-xcenter);
xmin = abs(xmin-xcenter);
minmax = xmin < xmax;

if (minmax)
  tmp = xmax;
  xmax = xmin;
  xmin = tmp;
end

% search for best logarithmic spacing
k = 0;
while k<=N
  q = ratio.^(-1/(N+k-1));
  d = xmax/(geometric_series(1,q,N+k-1) + s0);
  if d*(s0 + geometric_series(1,q,N+k-1)) <= xmin + eps
    break;
  end
  k = k+1;
end

x = d*(s0 + geometric_series(1,q,0:N+k-1));
if even
  x = [-x(N-k:-1:1), x];
else
  x = [-x(N-k:-1:1), 0, x];
end
if (minmax)
  x = fliplr(-x);
end
x = xcenter + x.';
end

function s = geometric_series(a0,q,N)
isargscalar(a0);
isargpositivescalar(q);
if (q == 1)
  s = a0.*(N+1);
else
  s = (1-q.^(N+1))./(1-q);
end
end
