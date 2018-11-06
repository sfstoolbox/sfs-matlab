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

switch area
  case 'position'
    conf = varargin{1};
    % the position x is a circular area with radius 0
    minmax_kGt_fun = @(x0p,xp) minmax_kt_circle(x0p,xp,0);  % handle
  case {'circle', 'circular'}
    Rl = varargin{1};  % radius of circular area
    conf = varargin{2};
    minmax_kGt_fun = @(x0p,xp) minmax_kt_circle(x0p,xp,Rl);  % handle
  case 'ellipse'
    to_be_implemented(mfilename);
  otherwise
    error('%s: unknown type (%s) for listening area', upper(mfilename), area);
end

N0 = conf.secondary_sources.number;
switch conf.secondary_sources.geometry
  case {'circle', 'circular'}
    R = conf.secondary_sources.size / 2;
    deltax = 2*pi*R/N0;
  case 'linear'
    L = conf.secondary_sources.size;
    deltax = L/(N0-1);
  otherwise
    error('%s: unsupported geometry (%s) for secondary sources', ...
      upper(mfilename), conf.secondary_sources.geometry);
end

% densely sampled SSD for calculating the aliasing frequency
conf.secondary_sources.number = 2048;
x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,xs,src);  % a_S(x0) >= 0
x0(:,7) = deltax;  % correct sampling

% local wavenumber of virtual sound field at x0
kSx0 = local_wavenumber_vector(x0(:,1:3), xs, src);

switch method
  case {'wfs', 'WFS'}
    % infinite control area
    minmax_kSt_fun = @(x0p) minmax_kt_circle(x0p,[0,0,0],Inf);  % handle
    f = aliasing_extended_control(x0, kSx0, x, minmax_kGt_fun, ...
      minmax_kSt_fun, conf);
  case {'nfchoa', 'HOA'}
    % circular control area at array center with Rl = M/k
    M = nfchoa_order(N0,conf);
    xc = conf.secondary_sources.center;
    f = aliasing_extended_modal(x0, kSx0, x, minmax_kGt_fun, xc, M, conf);
  case {'localwfs-sbl', 'LWFS-SBL'}
    % circular control area at xref with Rl = M/k
    xc = conf.xref;
    M = conf.localwfs_sbl.order;
    f = aliasing_extended_modal(x0, kSx0, x, minmax_kGt_fun, xc, M, conf);
  case {'localwfs-vss', 'LWFS-VSS'}
    switch conf.localwfs_vss.geometry
      case {'circle', 'circular'}
        xc = conf.localwfs_vss.center;
        Rc = conf.localwfs_vss.size / 2;
        minmax_kSt_fun = @(x0p) minmax_kt_circle(x0p,xc,Rc);  % handle
      otherwise
        error('%s: unsupported geometry (%s) for LWFS-VSS', ...
          upper(mfilename), conf.localwfs_vss.geometry);
    end
    f = aliasing_extended_control(x0, kSx0, x, minmax_kGt_fun, ...
      minmax_kSt_fun, conf);
  otherwise
    error('%s: unsupported SFS method (%s)', ...
      upper(mfilename), method)
end

