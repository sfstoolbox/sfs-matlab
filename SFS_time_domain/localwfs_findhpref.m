function [hpreflow,hprefhigh] = localwfs_findhpref(X,phi,xs,src,conf)
%LOCALWFS_FINDHPREF finds the frequency limits for the WFS pre-equalization
%filter
%
%   Usage: [hpreflow, hprefhigh] = localwfs_findhpref(X,phi,xs,src,conf)
%
%   Input parameters:
%       X       - listener position / m
%       phi     - listener direction [head orientation] / rad
%                 0 means the head is oriented towards the x-axis.
%       xs      - virtual source position / m
%       src     - source type: ...
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       hpreflow    - lower frequency limit of preequalization filter / Hz
%       hprefhigh   - higher frequency limit of preequalization filter / Hz
%
% See also: wfs_preequalization, wfs_fir_prefilter, wfs_iir_prefilter

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
nargmin = 5;
nargmax = 5;
narginchk(nargmin,nargmax);
if conf.debug
    isargposition(X);
    isargxs(xs);
    isargscalar(phi);
    isargchar(src);
    isargstruct(conf);
end


%% ===== Configuration ==================================================
conf.plot.useplot = false;    % disable plotting in spectrum_from_signal()
conf.ir.usehcomp = false;
conf.wfs.usehpre = false;     % no prefilter
conf.localwfs_vss.wfs.usehpre = false;  % no prefilter


%% ===== Variables ======================================================
N = conf.N;
irs = dummy_irs(round(N/2),conf);    % impulse responses
fs = conf.fs;                        % sampling rate
dimension = conf.dimension;          % dimensionality


%% ===== Computation ====================================================
% Compute impulse response/amplitude spectrum without prefilter
ir = ir_localwfs(X,phi,xs,src,irs,conf);
[H,~,f] = spectrum_from_signal(ir(:,1),conf);

H = H./H(1);  % normalize amplitude spectrum with H(f=0Hz)

% Model of local WFS spectrum without prefilter:
%   ^
% 1_| ______flow
%   |       \
%   |        \
%   |         \
%   |          fhigh
%   -------------------------> f

if strcmp('2.5D',dimension)  
    % Expected slope: 6dB per frequency-doubling
    % Find 6dB cut-off frequency
    flowidx = find(H <= 1/2, 1, 'first');
    hpreflow = f(flowidx)/2;
elseif strcmp('3D',dimension) || strcmp('2D',dimension)
    % Expected slope: 12dB per frequency-doubling
    % Find 12dB cut-off frequency 
    flowidx = find(H <= 1/4, 1, 'first');
    hpreflow = f(flowidx)/4;
else
    error('%s: %s is not a valid conf.dimension entry',upper(mfilename));
end

% Approximated slope beginning at hpreflow
Hslope = hpreflow./f(flowidx:end);
% Mean of H(f) evaluated from f to fs/2
Hmean = cumsum(H(end:-1:flowidx));   % cumulative sum
Hmean = Hmean./(1:length(Hmean)).';  % cumulative mean
Hmean = fliplr(Hmean);
% fhighidx is the frequency where both functions intersect the first time
fhighidx = flowidx - 1 + find( Hslope <= Hmean, 1, 'first');

if isempty(fhighidx)
  hprefhigh = fs/2;
else
  hprefhigh = f(fhighidx);
end 
