function xv = virtual_source_positions(x0,xs,src,conf)
%SECONDARY_SOURCE_POSITIONS Generates the positions and directions of the
%   virtual sources

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
isargsecondarysource(x0);
isargchar(src);
if nargin<nargmax
  conf = SFS_config;
else
  isargstruct(conf);
end

%% ===== Configuration ===================================================
virtualconf = conf;
virtualconf.secondary_sources.size = conf.virtual_secondary_sources.size;
virtualconf.secondary_sources.center = conf.virtual_secondary_sources.center;
virtualconf.secondary_sources.geometry = conf.virtual_secondary_sources.geometry;
virtualconf.secondary_sources.number = conf.virtual_secondary_sources.number;

geometry = conf.virtual_secondary_sources.geometry;
nls = virtualconf.secondary_sources.number;

%% ===== Main ============================================================

if strcmp('circle',geometry) || strcmp('circular',geometry)
  Rlocal = virtualconf.secondary_sources.size/2;  % radius of local area
  Xlocal = virtualconf.secondary_sources.center;  % center of local area

  % determine vector poiting towards source
  if strcmp('pw',src)
    ns = bsxfun(@rdivide,-xs,vector_norm(xs,2));
  else
    ns = bsxfun(@rdivide,xs-Xlocal,vector_norm(xs-Xlocal,2));
  end
  phis = atan2(ns(2),ns(1));  % azimuth angle of ns

  % == valid arc for virtual secondary sources (based on sec source position) ==
  delta_max = 0;
  delta_min = 0;
  % for each secondary source
  for idx=1:size(x0,1)
    xc0 = x0(idx,1:3) - Xlocal;  % vector from secondary source to local area
    Rc0 = vector_norm(xc0,2);  % distance from secondary source to local area
    nc0 = xc0./Rc0;  % normal vector from secondary source to local area
    % 1/2 opening angle of cone spanned by local area and secondary source
    phid = acos(Rlocal./Rc0);

    phiso = asin(ns(1)*nc0(2) - ns(2)*nc0(1));  % angle between ns and nc0
    delta_max = max(delta_max, phiso + phid);
    delta_min = min(delta_min, phiso - phid);
  end

  % == further constrain arc by position of virtual source ==
  if strcmp('pw',src)
    phid = pi/2;
  else
    % 1/2 opening angle of cone spanned by local area and virtual source
    phid = acos(Rlocal./vector_norm(xs-Xlocal,2));
  end

  delta_max = min(delta_max, phid);
  delta_min = max(delta_min, -phid);
  delta_offset = (delta_max - delta_min)/(2*nls);

  % === equi-distant sampling on valid arc ===
  % Azimuth angles
  phi = phis + linspace(delta_min + delta_offset,delta_max -delta_offset, nls)';
  % Elevation angles
  theta = zeros(nls,1);
  % Positions of the secondary sources
  [cx,cy,cz] = sph2cart(phi,theta,Rlocal);
  xv(:,1:3) = [cx,cy,cz] + repmat(Xlocal,nls,1);
  % Direction of the secondary sources
  xv(:,4:6) = direction_vector(xv(:,1:3),repmat(Xlocal,nls,1).*ones(nls,3));
  % equal weights for all sources
  xv(:,7) = ones(nls,1);
else
  xv = secondary_source_positions(virtualconf);
  xv = secondary_source_selection(xv,xs,src);
end

xv = secondary_source_tapering(xv,virtualconf);
