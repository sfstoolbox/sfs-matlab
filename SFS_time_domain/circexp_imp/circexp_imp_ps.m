function [pm,delay_offset] = circexp_imp_ps(xs,Nce,xq,fhp,conf)
%CIRCEXP_IMP_PS calculates the circular basis expansion of a point source
%
%   Usage: [pm,delay_offset] = circexp_imp_ps(xs,Nce,xq,[fhp],conf)
%
%   Input parameters:
%       xs      - position of point source / m [1 x 3]
%       Nce     - maximum order of circular basis expansion
%       xq      - optional expansion center / m [1 x 3]
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       pm            - regular circular expansion coefficients in time domain
%                       for m = 0:Nce, [conf.N x Nce+1]
%       delay_offset  - additional added delay, so you can correct it
%
%   CIRCEXP_IMP_PS(xs,Nce,xq,fhp,conf) returns the circular basis expansion of
%   a point source.
%
%   See also: circexp_imp_pw

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Developers                             *
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


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
isargcoord(xq,xs);
isargpositivescalar(Nce);
if nargin == nargmin
  conf = fhp;
else
  isargpositivescalar(fhp);
end


%% ===== Configuration ==================================================
N = conf.N;
c = conf.c;
fs = conf.fs;


%% ===== Computation =====================================================

xs = xs - xq; % shift coordinates
[phis, rs] = cart2pol(xs(1),xs(2));

%-------------------------------------------------------------------------------
% Implementation of
%
%          j^(|m|-m)
% p_m(t) = --------- IFT[ -jk h_|m|(k rs) ] e^(-im phis) e^(im phi0)
%             4pi
%
% with IFT being the inverse fourier transform of the spherical hankel function.
% The IFT is realised by an IIR-implementation using the Laplace domain (s)
% respresentation of the spherical hankel function:
%                                     _____
%                      exp(-j s tau)   | |  (s - z_n/tau)
% h_|m|(s tau) = j^|m| -------------   | |  -------------
%                        -j s tau       n   (s - p_n/tau)
%
% z_i and p_i are the zeros and poles of the hankel function. tau = r_s/c.
%-------------------------------------------------------------------------------

% === Linkwitz-Riley (LR) filter for stabilisation ===
zlr = [];
plr = [];
klr = 1;
if Nce > 0
  if nargin == nargmin
    warning('without highpass filtering the implementation is not stable');
  else
    [zlr, plr, klr] = linkwitz_riley(ceil(Nce/2)*2, fhp/fs*2, 'high');
  end
end

% Compute impulse responses for each mode m
pulse = dirac_imp();
pm = [repmat(pulse,[1 Nce+1]); zeros(N-length(pulse),Nce+1)];
% Negative m can be inferred from symmetry relations
for m=0:Nce
    % === IIR-Implementation of Spherical Hankel function ===
    if m==0  % not supported by octave
        zh = []; ph = []; kh=1;
    else 
        [zh, ph] = sphbesselh_zeros(m);
        % bilinear transform to z-domain
        if isoctave
            [zh, ph, kh] = bilinear(zh*c/rs, ph*c/rs, 1, 1/fs);
        else
            [zh, ph, kh] = bilinear(zh*c/rs, ph*c/rs, 1, fs);
        end
    end   
    % === Apply Hankel + LR Filter to current mode ===
    % zeros remaining after compensating the poles of hankel function (ph)
    zlr_comp = ones(length(zlr)-length(ph),1);
    % generate second-order-sections
    [sos, g] = zp2sos([zlr_comp; zh], plr, klr*kh, 'down', 'none');
    % filtering
    pm(:,m+1) = sosfilt(sos, pm(:,m+1)).*g.*(1j).^m.*exp(-1j*m*phis);
end
pm = pm./(4*pi*rs);

delay_offset = -rs/c;
