function [d,delay_offset] = driving_function_imp_localwfs_sbl_ps(x0,xs,conf)

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

%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 3;
narginchk(nargmin,nargmax);

%% ===== Configuration ========================================================
N0 = size(x0,1);
Nfft = conf.N;
xref = conf.xref;
fs = conf.fs;
c = conf.c;
% Ambisonics order
if isempty(conf.localsfs.sbl.order)
    Nce = nfchoa_order(N0,conf);
else
    Nce = conf.localsfs.sbl.order;
end
% Crossover frequency
if isempty(conf.localsfs.sbl.fc)
  fc = aliasing_frequency(conf);
else
  fc = conf.localsfs.sbl.fc;  
end
% Resolution of plane wave decomposition
if isempty(conf.localsfs.sbl.Npw)
  Npw = 2*ceil(2*pi*0.9*fs/conf.c*conf.secondary_sources.size/2);
else
  Npw = conf.localsfs.sbl.Npw;
end
driving_functions = conf.driving_functions;

conf.t0 = 'source';  % needed to time-align lf and hf part of driving function
wfsconf = conf;
wfsconf.wfs = conf.localsfs.wfs;
wfsconf.driving_functions = conf.localsfs.sbl.driving_functions;

%% ===== Computation ==========================================================

% === Local WFS for high frequencies ===
switch driving_functions
case {'default', 'hp-filtered'}
  % regular circular expansion of point source (highpass implicitly)
  [pm, delay_circexp] = circexp_imp_ps(xs, Nce, xref, fc, conf);
  % coefficient for lowpass filtering of WFS driving function
  [zlp, plp, klp] = linkwitz_riley(ceil(Nce/2)*2, fc/fs*2, 'low');  % lr-filter
case 'hf-approximation'
  % regular circular expansion of point source (high-frequency approximation)
  rs = norm(xs - xref);
  ns = -(xs-xref)./rs;
  %
  pm = circexp_imp_pw(ns, Nce, xref, conf);
  pm = pm/(4*pi*rs);
  %
  delay_circexp = -rs/c;
  % highpass filtering
  [zhp, php, khp] = butter(1, fc/fs*2, 'high');
  [sos, g] = zp2sos(zhp, php, khp, 'down', 'none');  % generate sos
  pm = sosfilt(sos, pm, 2)*g;
  % coefficient for lowpass filtering of WFS driving function
  [zlp, plp, klp] = butter(1, fc/fs*2, 'low');
end
% plane wave decomposition
ppwd = pwd_imp_circexp(pm, Npw);
% driving signal
[d_lwfs, delay_lwfs] = driving_function_imp_wfs_pwd(x0, ppwd, xref, wfsconf);

% === WFS for low frequencies ===
% 
conf.driving_functions = 'reference_point';
% secondary source selection
[~, xdx] = secondary_source_selection(x0, xs, 'ps');
x0(xdx,:) = secondary_source_tapering(x0(xdx,:) , conf);
% driving function
d_lp = zeros(Nfft,N0);
[d_lp(:,xdx), ~, ~, delay_lp] = driving_function_imp_wfs(x0(xdx,:), xs, 'ps', conf);
% lowpass filtering
[sos, g] = zp2sos(zlp, plp, klp, 'down', 'none');  % generate sos
d_lp = sosfilt(sos, d_lp, 1)*g;
  
% === Crossover ===
% get delay of delayline
[~, delay_delayline] = delayline(1, 0, 0, conf);
% delay to compensate between lf-part and hf-part
delay_comp = (delay_lp-delay_lwfs-delay_circexp)*conf.fs - delay_delayline;
% final driving signal
d = d_lp + delayline(d_lwfs, delay_comp, 1, conf);
% final delay
delay_offset = delay_lp;
