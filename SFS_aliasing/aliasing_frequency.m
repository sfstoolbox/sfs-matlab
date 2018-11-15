function f = aliasing_frequency(method, xs, src, area, x, varargin)
%ALIASING_FREQUENCY Summary of this function goes here
%   Detailed explanation goes here

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2018 SFS Toolbox Developers                             *
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

%% ===== Checking input parameters =======================================
nargmin = 6;
nargmax = 7;
narginchk(nargmin,nargmax);
isargchar(method, area, src);
isargxs(xs);
ismatrix(x);

switch area
  case 'position'
    % the position x is a circular area with radius 0
    minmax_kGt_fun = @(x0p,xp) minmax_kt_circle(x0p,xp,0);  % handle
    conf = varargin{1};
  case 'line'
    to_be_implemented(mfilename);
  case {'circle', 'circular'}
    Rl = varargin{1};  % radius of circular area
    isargpositivescalar(Rl);
    minmax_kGt_fun = @(x0p,xp) minmax_kt_circle(x0p,xp,Rl);  % handle
    conf = varargin{2};
  case 'ellipse'
    to_be_implemented(mfilename);
  otherwise
    error('%s: unknown type (%s) for listening area', upper(mfilename), area);
end
isargstruct(conf);   

%% ===== Computation ====================================================

% prepare virtual secondary sources
N0 = conf.secondary_sources.number;
switch conf.secondary_sources.geometry
  case {'circle', 'circular'}
    R = conf.secondary_sources.size / 2;
    deltax0 = 2*pi*R/N0;
  case 'linear'
    L = conf.secondary_sources.size;
    deltax0 = L/(N0-1);
  otherwise
    error('%s: unsupported geometry (%s) for secondary sources', ...
      upper(mfilename), conf.secondary_sources.geometry);
end

% densely sampled SSD for calculating the aliasing frequency
conf.secondary_sources.number = 256;
x0 = secondary_source_positions(conf);
x0(:,7) = deltax0;  % correct sampling

% normalised local wavenumber vector of virtual sound field kS(x0)
kSx0 = local_wavenumber_vector(x0(:,1:3), xs, src);

switch method
  case {'wfs', 'WFS'}
    M = Inf;  % WFS does not have any bandwidth limitation
    xc = [0,0,0];  % no effect since M = inf
    Npw = Inf;  % WFS does not involve a discretised Plane Wave Decomposition
    f = aliasing_extended_sbl(x0,kSx0,x,minmax_kGt_fun,xc,M,Npw,conf);
  case {'nfchoa', 'HOA'}
    % circular control area at array center with Rl = M/k
    M = nfchoa_order(N0,conf);  % 
    xc = conf.secondary_sources.center;  % expansion around the array center
    Npw = Inf;  % HOA does not involve a discretised Plane Wave Decomposition
    f = aliasing_extended_sbl(x0,kSx0,x,minmax_kGt_fun,xc,M,Npw,conf);
  case {'localwfs-sbl', 'LWFS-SBL'}
    % circular control area at xref with Rl = M/k
    M = conf.localwfs_sbl.order;
    xc = conf.xref;
    Npw = conf.localwfs_sbl.Npw;
    f = aliasing_extended_sbl(x0,kSx0,x,minmax_kGt_fun,xc,M,Npw,conf);
  case {'localwfs-vss', 'LWFS-VSS'}
    
    % prepare virtual secondary sources
    Nv = conf.localwfs_vss.number;
    switch conf.localwfs_vss.geometry
      case {'circle', 'circular'}
        Rv = conf.localwfs_vss.size / 2;
        deltaxv = 2*pi*Rv/Nv;
      case 'linear'
        Lv = conf.localwfs_vss.size;
        deltaxv = Lv/(Nv-1);
      otherwise
        error('%s: unsupported geometry (%s) for virtual secondary sources', ...
          upper(mfilename), conf.localwfs_vss.geometry);
    end
    
    % densely sampled virtual SSD for calculating the aliasing frequency
    virtualconf = conf;
    virtualconf.secondary_sources.size = conf.localwfs_vss.size;
    virtualconf.secondary_sources.center = conf.localwfs_vss.center;
    virtualconf.secondary_sources.geometry = conf.localwfs_vss.geometry;
    virtualconf.secondary_sources.number = 256;
    
    xv = secondary_source_positions(virtualconf);
    xv(:,7) = deltaxv;
    
    % local wavenumber of virtual sound field at xv
    kSxv = local_wavenumber_vector(xv(:,1:3), xs, src);

    f = aliasing_extended_vss(x0,xv,kSxv,x,minmax_kGt_fun,conf);
  otherwise
    error('%s: unsupported SFS method (%s)', ...
      upper(mfilename), method)
end

