function f = aliasing_extended_sbl(x0, kS, x, minmax_kGt_fun, xc, M, Npw, conf)
%ALIASING_EXTENDED_SBL aliasing frequency for an extended listening area for
%an circular control area at xc with R=M/k where synthesis is focused on.
%
%   Usage: f = aliasing_extended_control(x0, kSx0, x, minmax_kGt_fun, minmax_kSt_fun, conf)
%
%   Input options:
%       x0              - position, direction, and sampling distance of
%                         secondary sources [N0x7] / m
%       kSx0            - normalised local wavenumber vector of virtual sound
%                         field [N0x3]
%       x               - position for which aliasing frequency is calculated
%                         [Nx3]
%       minmax_kGt_fun  - function handle to determine the extremal value of
%                         the tangential component of k_G(x-x0)
%                         [kGtmin, kGtmax] = minmax_kGt_fun(x0)
%       xc              - center of circular control area
%       M               - modal order which defines the radius R=M/k
%       Npw             - 
%       conf            - configuration struct (see SFS_config)
%
%   Output parameters:
%       f   - aliasing frequency [Nx1]
%

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

% shift coordinates to expansion centre xc
x = bsxfun(@minus, x, xc);
x0(:,1:3) = bsxfun(@minus, x0(:,1:3), xc);

phin0 = cart2pol(x0(:,4),x0(:,5));  % azimuth angle of normal vector n0(x0)
deltax0 = abs(x0(:,7));  % sampling distance delta_x0(x0)
phik = cart2pol(kS(:,1),kS(:,2));  % azimuth angle of kS(x0)

select = cos(phin0 - phik) >= 0;  % secondary source selection
gdx = find(select);

% quantities for selected secondary sources
x0select = x0(select,:);
[phix0select, rx0select] = cart2pol(x0select(:,1),x0select(:,2));
deltax0select = deltax0(select);  % sampling distance
kSselect = kS(select,:);
phikselect = phik(select);
kStselect = sin(phin0(select) - phikselect); % tangential component of kS(x0)

% wavelength at which secondary source x0 turn inactive due to spatial
% bandwidth limitation
lambdaMselect = 2.*pi/M.*rx0select.*abs(sin(phix0select - phikselect));

f = inf(size(x,1), 1);
if Npw < M
  f(:) = 0;
  return
end

eps = 0.01;

for xdx = 1:size(x,1);
  
  vectorxx0 = bsxfun(@minus, x(xdx,:), x0(:,1:3));  % vector x - x0
  [phixx0,rxx0] = cart2pol(vectorxx0(:,1),vectorxx0(:,2));  % polar coordinates
  % mininum and maximum values of k_Gt(x - x_0) 
  % (tangential component of k_G(x-x0))
  [kGtmin, kGtmax] = minmax_kGt_fun(x0,x(xdx,:));
  kGtminselect = kGtmin(gdx);
  kGtmaxselect = kGtmax(gdx);

  lambda = 0;
  
  % aliasing wavelength for eta=1 and zeta=0 
  lambda10 = deltax0select.* ...
    max(abs(kStselect-kGtminselect),abs(kStselect-kGtmaxselect));
  % select all wavelengths, that above the limit of the spatial bandwidth 
  % limitation
  select = lambda10 >= lambdaMselect;
  if any(select)
    lambda = max(lambda, max(lambda10(select)));
  end
  
  if isinf(Npw)
    f(xdx) = conf.c./lambda;
    continue;
  end  
  
  % aliasing wavelength for eta=0 and zeta=1
  lambda01 = 2.*pi./Npw.*rxx0(gdx).*abs(sin(phixx0(gdx)-phikselect));
  % select all wavelengths, that above the limit of the spatial bandwidth 
  % limitation
  select = lambda01 >= lambdaMselect;
  if any(select)
    lambda = max(lambda, max(lambda01(select)));
  end
  
  % aliasing wavelength for eta=1 and zeta=1
  for sdx=1:size(x0select,1);  % interate over all active x0'
    
    % secondary source selection for ALL x0 with a single kS(x0')
    active = x0(:,4:6)*kSselect(sdx,:).' >= 0;
    
    % vector x0 - x0'
    vectorx0x0p = bsxfun(@minus, x0(active,1:3), x0select(sdx,1:3));
    % polar angle and radius of vector (x0-x0')
    [phix0x0p,rx0x0p] = cart2pol(vectorx0x0p(:,1),vectorx0x0p(:,2));
    
    % range for aliasing wavelength at x0' caused by Discrete SSD (eta = 1)
    lambda10x0p1 = ...
      abs(sin(phin0(active) - phikselect(sdx))-kGtmin(active)) ...
      .*deltax0(active);
    lambda10x0p2 = ...
      abs(sin(phin0(active) - phikselect(sdx))-kGtmax(active)) ...
      .*deltax0(active);
    lambda10x0pmin = min(lambda10x0p1,lambda10x0p2);
    lambda10x0pmax = max(lambda10x0p1,lambda10x0p2);
    
    % aliasing wavelength at x0' caused by Plane Wave Decomposition (zeta = 1)
    lambda01x0p = 2.*pi/Npw.*rx0x0p.*abs(sin(phix0x0p - phikselect(sdx)));
    
    % select x0' where PWD wavelength is in the range of SSD wavelength
    select = lambda01x0p < (1+eps).*lambda10x0pmax & ...
        lambda10x0pmin < (1+eps).*lambda01x0p;
    if any(select)
      lambda11 = lambda01x0p(select);
      % select all wavelengths, that above the limit of the spatial bandwidth 
      % limitation
      select = lambda11 >= lambdaMselect(sdx);
      if any(select)
        lambda = max(lambda, max(lambda11(select)));
      end
    end
  end  

  f(xdx) = conf.c./lambda;
end
