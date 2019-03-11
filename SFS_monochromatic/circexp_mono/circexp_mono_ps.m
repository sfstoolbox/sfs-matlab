function Pm = circexp_mono_ps(xs,Nce,f,xq,fhp,conf)
%CIRCEXP_MONO_PS circular basis expansion of a mono-frequent point source
%
%   Usage: Pm = circexp_mono_ps(xs,Nce,f,xq,[fhp],conf)
%
%   Input parameters:
%       xs      - position of point source / m [1 x 3]
%       Nce     - maximum order of circular basis expansion
%       f       - frequency of the monochromatic source / Hz
%       xq      - expansion center / m [1 x 3]
%       fhp     - cut-off frequency of highpass for regularisation / Hz
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       Pm      - regular circular expansion coefficients
%                 for m = 0:Nce, [1 x Nce+1]
%
%   See also: circexp_mono_pw

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2019 SFS Toolbox Developers                             *
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
% https://sfs.readthedocs.io                            sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 6;
narginchk(nargmin,nargmax);
isargcoord(xq,xs);
isargpositivescalar(Nce,f);
if nargin == nargmin
  conf = fhp;
else
  isargpositivescalar(fhp);
end

%% ===== Configuration ==================================================
c = conf.c;


%% ===== Computation =====================================================
xs = xs - xq;  % shift coordinates
[phis, rs] = cart2pol(xs(1),xs(2));
tau = rs/c;
omega = 2*pi*f;
omega_vec = [-omega.^2; 1i*omega; 1];  % vector for evaluation of SOS

%-------------------------------------------------------------------------------
% Implementation of
%
%                       i^(|m|-m)
% P_m = -ik h_|m|(k rs) --------- e^(-im phis)
%                          4pi
%
% The Laplace domain (s) respresentation of the spherical Hankel function:
%                                     _____
%                      exp(-s tau)     | |  (s - z_n/tau)
% h_|m|(s tau) = i^|m| -------------   | |  -------------
%                        -i s tau       n   (s - p_n/tau)
%
% z_n and p_n are the zeros and poles of the Hankel function. tau = r_s/c.
%-------------------------------------------------------------------------------

% === Linkwitz-Riley (LR) filter for stabilisation ===
zlr = [];
plr = [];
klr = 1;
if Nce > 0
  if nargin == nargmin
    warning('without highpass filtering the implementation is not stable');
  else
    % zero-pole-gain in s-domain of LR-Filter
    Nlr = ceil(Nce/2)*2;
    [zlr, plr, klr] = linkwitz_riley(Nlr, 2*pi*fhp, 'high', 's');
  end
end

% Compute values for each mode m
Pm = zeros(1,2*Nce+1);
for m=0:Nce  % Negative m can be inferred from symmetry relations
    % === zero-pole-gain in s-domain of Spherical Hankel function ===
    kh=1;
    if m==0  % not supported by octave
        zh = []; ph = [];
    else 
        [zh, ph] = sphbesselh_zeros(m);
        zh = zh./tau;
        ph = ph./tau;        
    end
    % Zeros/poles remaining after compensating the poles of Hankel function
    zlr_comp = zeros(length(zlr)-length(ph),1);  % empty if negative
    ph_comp = zeros(length(ph)-length(zlr),1);  % empty if negative
    % Generate second-order-sections
    [sos, g] = zp2sos([zlr_comp; zh], [ph_comp; plr], klr*kh, 'down', 'none');
    % Compute value of sos at s = iw
    Pm(m+Nce+1) = g.*prod( (sos(:,1:3)*omega_vec)./(sos(:,4:6)*omega_vec),1);
end

Pm(1:Nce) = Pm(2*Nce+1:-1:Nce+2);  % symmetry
Pm = Pm.*exp(+1i*(-Nce:Nce)*(pi/2-phis));  % apply terms only depending on m
Pm = Pm./(4*pi*rs).*exp(-1i*omega*tau);  % scaling and delay
