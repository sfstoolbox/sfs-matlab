function [d,delay_offset] = driving_function_imp_wfs_vss(x0,xv,srcv,dv,conf)
%DRIVING_FUNCTION_IMP_WFS_VSS returns the driving signal d for a given set of
%virtual secondary sources and their corresponding driving signals
%
%   Usage: d = driving_function_imp_wfs_vss(x0,dv,xv,srcv,conf)
%
%   Input parameters:
%       x0          - position, direction, and weights of the real secondary
%                     sources / m [N0x7]
%       xv          - position, direction, and weights of the virtual secondary
%                     sources / m [Nvx7]
%       srcv        - type of virtual secondary sources
%       dv          - driving signals of virtual secondary sources [NxNv]
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       d             - driving function signal [NxN0]
%       delay_offset  - additional added delay, so you can correct it
%
%   References:
%       S. Spors (2010) - "Local Sound Field Synthesis by Virtual Secondary
%                          Sources", 40th AES
%
%   See also: driving_function_imp_localwfs, driving_function_mono_wfs_vss

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
isargmatrix(dv);
isargsecondarysource(x0,xv);
isargchar(srcv);
isargstruct(conf);

%% ===== Computation ====================================================
% Initialize
Nv = size(xv,1);
N0 = size(x0,1);

tau0 = inf(Nv,N0);
w0 = zeros(Nv,N0);

% Calculate pre-equalization filter if required hpre(t)
[hpre,delay_offset] = wfs_preequalization(dirac_imp(),conf);

% apply hpre to pwd if its more efficient (less virtual sources than 
% loudspeakers)
if Nv <= N0
    dv = convolution(dv,hpre);
end
N = size(dv,1);

% secondary source selection and driving function to synthesise a single virtual
% secondary source
switch srcv
case 'fs'
    ssd_select = @(X0,XS) secondary_source_selection(X0,XS(1:6)','fs');
    driv = @(X0,XS) driving_function_imp_wfs_fs(X0(:,1:3),X0(:,4:6),XS,conf);
case 'pw'
    ssd_select = @(X0,XS) secondary_source_selection(X0,XS(1:3)','pw');
    driv = @(X0,XS) driving_function_imp_wfs_pw(X0(:,1:3),X0(:,4:6),XS,conf);
end

idx = 1;
for xvi = xv'
    % select active source for single virtual secondary source
    [x0s, xdx] = ssd_select(x0,xvi);
    if ~isempty(x0s) && xvi(7) > 0
        % Focused source position
        xs = repmat(xvi(1:3)',[size(x0s,1) 1]);
        % Delay and weights for single virtual secondary source
        [tau0(idx,xdx),w0(idx,xdx)] = driv(x0s,xs);
        % Optional tapering
        wtap = tapering_window(x0s, conf);
        % Apply secondary sources' tapering and possibly virtual secondary
        % sources' tapering to weighting matrix
        w0(idx,xdx) = w0(idx,xdx).*wtap.'.*xvi(7);
    end
    idx = idx + 1;
end

% Remove delay offset, in order to begin always at t=0 with the first wave 
% front at any secondary source
delay_offset = delay_offset - min(tau0(w0~=0));
tau0 = tau0 - min(tau0(w0~=0)); 

% Compose impulse responses
d = zeros(N,N0);
for ndx=1:Nv
    % select only secondary source with non-zero weights
    xsel = w0(ndx,:) ~= 0;
    
    [tmp, delayline_delay] = ...
        delayline(dv(:,ndx),tau0(ndx,xsel),w0(ndx,xsel),conf);
    d(:,xsel) = d(:,xsel) + tmp;
end
% add delay of delayline
delay_offset = delay_offset + delayline_delay;

% apply hpre to driving signals if its more efficient
if Nv > N0
    d = convolution(d,hpre);
end
d = d(1:conf.N, :);
