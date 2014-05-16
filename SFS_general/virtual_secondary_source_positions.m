function xv = virtual_secondary_source_positions(x0,xs,src,conf)
%VIRTUAL_SECONDARY_SOURCE_POSITIONS Generates the positions and directions of the
%   virtual sources
%
%   Usage: xv = virtual_secondary_source_positions(x0,xs,src,conf)
%
%   Input options:
%       x0          - positions, directions and weights of real secondary
%                     sources [nx7]
%       xs          - position and for focused sources also direction of the
%                     desired source model / m [1x3] or [1x6]
%       src         - source type of the virtual source
%                       'pw' - plane wave (xs is the direction of the
%                              plane wave in this case)
%                       'ps' - point source
%                       'fs' - focused source
%       conf        - optional configuration struct (see SFS_config)
%
%   Output options:
%       xv          - virtual secondary source positions, directions and
%                     weights / m
%
%   VIRTUAL_SECONDARY_SOURCE_POSITIONS(x0,xs,src,conf) generates the positions
%   and directions xv of virtual secondary sources for a given geometry
%   (conf.localsfs.geometry) and array size (conf.localsfs.size).
%   The direction of the virtual sources is given as their unit vectors
%   pointing in the given direction. The algorithm determines the (optimal)
%   positioning of the virtual secondary sources by taking the positions of
%   the virtual source and the real sources into account.

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
virtualconf.secondary_sources.size = conf.localsfs.size;
virtualconf.secondary_sources.center = conf.localsfs.center;
virtualconf.secondary_sources.geometry = conf.localsfs.geometry;
virtualconf.secondary_sources.number = conf.localsfs.vss.number;

geometry = conf.localsfs.geometry;
sampling = conf.localsfs.vss.sampling;
nls = conf.localsfs.vss.number;
ignore_secondary_sources = conf.localsfs.vss.ignoress;

%% ===== Main ============================================================

if strcmp('circle',geometry) || strcmp('circular',geometry)
  % =====================================================================
  % virtual circular Array
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

  % CONSTRAINT 1 ========================================================
  % valid arc by position of virtual source
  if strcmp('pw',src)
    phid = pi/2;
  else
    % 1/2 opening angle of cone spanned by local area and virtual source
    phid = acos(Rl./vector_norm(xs-xl,2));
  end

  % the positions of the 'real' secondary source can be taken into account
  % virtual secondary source position
  if (ignore_secondary_sources || isempty(x0))
    delta_max = phid;
    delta_min = -phid;
  else
    % CONSTRAINT 2 ========================================================
    % valid arc for virtual secondary sources (based on sec source positions)
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
    case 'linear'
      % === equi-angular sampling on valid arc ===
      phi = phis + linspace(delta_min + delta_offset,delta_max-delta_offset, nls).';
    case 'log'
      logratio = conf.localsfs.vss.logratio;

      % === logarithmic-angular sampling on valid arc ===
      if mod(nls,2) == 0
        N = nls./2;
        a0 = -0.5;
      else
        N = (nls-1)./2;
        a0 = 0.0;
      end

      phimax = max(abs(delta_min + delta_offset),abs(delta_max-delta_offset))
      phimin = min(abs(delta_min + delta_offset),abs(delta_max-delta_offset));

      k = 0;
      while k<=N
        q = logratio.^(-1/(N+k-1));
        d = phimax/(geometric_series(1,q,N+k-1) + a0);
        if d*(a0 + geometric_series(1,q,N+k-1)) <= phimin + eps
          break;
        end
        k = k+1;
      end

      phi = d*(a0 + geometric_series(1,q,0:N+k-1))
      if mod(nls,2) == 0
        phi = [-phi(N-k:-1:1), phi];
      else
        phi = [-phi(N-k:-1:1), 0, phi];
      end
      if (abs(delta_min + delta_offset) > abs(delta_max-delta_offset))
        phi = fliplr(-phi);
      end
      phi = phis + phi';
    case 'scalar'
      x = linspace(sin(delta_min + delta_offset),sin(delta_max - delta_offset),nls).';
      phi = phis + asin(x);
    case 'projective'
      x = linspace(-3,3,nls).';
      phi = phis + atan(x);
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
else
  % ===== traditional way =====

  % just select the virtual secondary sources based on the position and type
  % of the virtual source
  xv = secondary_source_positions(virtualconf);
  xv = secondary_source_selection(xv,xs,src);
end

% xv = secondary_source_tapering(xv,virtualconf);
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
