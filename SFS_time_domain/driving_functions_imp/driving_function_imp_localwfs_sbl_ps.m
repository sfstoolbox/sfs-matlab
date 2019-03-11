function [d,delay_offset] = driving_function_imp_localwfs_sbl_ps(x0,xs,conf)
%DRIVING_FUNCTION_IMP_LOCALWFS_SBL_PS driving signal for a point source using
%local WFS with spatial bandwidth limitation
%
%   Usage: [d,delay_offset] = driving_function_imp_localwfs_sbl_ps(x0,xs,conf)
%
%   Input parameters:
%       x0          - position and direction of the secondary source / m [nx7]
%       xs          - position of virtual point source / m [1x3]
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       d               - driving signals [mxn]
%       delay_offset    - additional added delay, so you can correct it
%
%   See also: sound_field_imp_localwfs_sbl, driving_function_imp_localwfs_sbl
%
%   References:
%       Winter, Hahn, Spors (2017), "Time-Domain Realisation of Model-Based
%       Rendering for 2.5D Local Wave Field Synthesis Using Spatial
%       Bandwidth-Limitation", 25th European Signal Processing Conference
%       (EUSIPCO), pp. 718-722, https://bit.ly/2yMjlOw 

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


%% ===== Checking of input  parameters ========================================
nargmin = 3;
nargmax = 3;
narginchk(nargmin,nargmax);


%% ===== Configuration ========================================================
N0 = size(x0,1);
Nfft = conf.N;
xref = conf.xref;
fs = conf.fs;
% Ambisonics order
if isempty(conf.localwfs_sbl.order)
    Nce = nfchoa_order(N0,conf);
else
    Nce = conf.localwfs_sbl.order;
end
% Crossover frequency
if isempty(conf.localwfs_sbl.fc)
    fc = aliasing_frequency(conf);
else
    fc = conf.localwfs_sbl.fc;
end
% Resolution of plane wave decomposition
if isempty(conf.localwfs_sbl.Npw)
    Npw = 2*ceil(2*pi*0.9*fs/conf.c*conf.secondary_sources.size/2);
else
    Npw = conf.localwfs_sbl.Npw;
end

%% ===== Variables ============================================================
Nlr = ceil(Nce/2)*2;  % order of Linkwitz-Riley Coefficients
Wlr = fc/fs*2;  % normalised cut-off frequency of Linkwitz-Riley


%% ===== Computation ==========================================================

% === Local WFS for high frequencies ===
% Regular circular expansion of point source (highpass implicitly)
% Winter et al. (2017), eq. (14)
[pm,delay_circexp] = circexp_imp_ps(xs,Nce,xref,fc,conf);
% Modal window
wm = modal_weighting(Nce,conf);
pm = bsxfun(@times,wm,pm);
% Plane wave decomposition, Winter et al. (2017), eq. (10)
ppwd = pwd_imp_circexp(pm,Npw);
% Driving signal, Winter et al. (2017), eq. (9)
[dlwfs,delay_lwfs] = driving_function_imp_wfs_pwd(x0,ppwd,xref,conf);

% === WFS for low frequencies ===
% Secondary source selection
[~,xdx] = secondary_source_selection(x0,xs,'ps');
x0(xdx,:) = secondary_source_tapering(x0(xdx,:),conf);
% Driving signal
dwfs = zeros(Nfft,N0);
[dwfs(:,xdx),~,~,delay_lp] = driving_function_imp_wfs(x0(xdx,:),xs,'ps',conf);
% Coefficients for lowpass filtering of WFS driving function
[zlp,plp,klp] = linkwitz_riley(Nlr,Wlr,'low');  % LR Lowpass filter
% Lowpass filtering
[sos,g] = zp2sos(zlp,plp,klp,'down','none');  % generate sos
dwfs = sosfilt(sos,dwfs)*g;
  
% === Crossover ===
% Winter et al. (2017), eq. (15)
% Get delay of delayline
[~,delay_delayline] = delayline(1,0,0,conf);
% Delay to compensate between lf-part and hf-part
delay_comp = delay_lp - (delay_lwfs+delay_circexp+delay_delayline);
% Combined driving signal
d = dwfs + delayline(dlwfs,delay_comp,1,conf);

% === Compensate Phase-Distortions by Inverse Allpass ===
% Winter et al. (2017), eq. (17)
% TODO: ensure that the time-reserved filter is not truncated
[zap,pap,kap] = linkwitz_riley(Nlr,Wlr,'all');  % LR Allpass filter
[sos,g] = zp2sos(zap,pap,kap,'down','none');  % generate sos
% (time-reversed) allpass filtering
d = sosfilt(sos,d(end:-1:1,:))*g;
d = d(end:-1:1,:);

% Final delay
delay_offset = delay_lp;
