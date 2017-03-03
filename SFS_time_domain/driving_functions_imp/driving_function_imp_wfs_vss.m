function d = driving_function_imp_wfs_vss(x0,xv,dv,conf)
%DRIVING_FUNCTION_IMP_WFS_VSS returns the driving signal d for a given set of
%virtual secondary sources and their corresponding driving signals
%
%   Usage: d = driving_function_imp_wfs_vss(x0,xv,dv,conf)
%
%   Input parameters:
%       x0          - position, direction, and weights of the real secondary
%                     sources / m [nx7]
%       xv          - position, direction, and weights of the virtual secondary
%                     sources / m [mx7]
%       dv          - driving signals of virtual secondary sources [Sxm]
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       d           - driving function signal [Sxn]
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
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargmatrix(dv);
isargsecondarysource(x0,xv);
isargstruct(conf);


%% ===== Configuration ==================================================
fs = conf.fs;
N = conf.N;


%% ===== Computation ====================================================
% Apply wfs preequalization filter on each driving signal of the vss'
dv = wfs_preequalization(dv, conf);

% Initialize
Nv = size(xv,1);
N0 = size(x0,1);

d = zeros(N,N0);
delay = inf(N0,Nv);
weight = zeros(N0,Nv);

idx = 1;
for xvi = xv'
    % Select active source for one focused source
    [x0s, xdx] = secondary_source_selection(x0,xvi(1:6)','fs');
    if ~isempty(x0s) && xvi(7) > 0
        % Focused source position
        xs = repmat(xvi(1:3)',[size(x0s,1) 1]);
        % Delay and weights for single focused source
        [delay(xdx,idx),weight(xdx,idx)] = ...
            driving_function_imp_wfs_fs(x0s(:,1:3),x0s(:,4:6),xs,conf);
        % Optional tapering
        x0s = secondary_source_tapering(x0s,conf);
        % Apply secondary sources' tapering and possibly virtual secondary
        % sources' tapering to weighting matrix
        weight(xdx,idx) = weight(xdx,idx).*x0s(:,7).*xvi(7);
    end
    idx = idx + 1;
end

% Remove delay offset, in order to begin always at t=0 with the first wave front
% at any secondary source
delay = delay - min(delay(:));

% Compose impulse responses
for idx=1:Nv
    xdx = weight(:,idx) ~= 0;
    if sum(xdx) > 0
        % Shift and weight prototype driving function
        pulse = repmat(dv(:,idx), 1, sum(xdx));
        d(:,xdx) = d(:,xdx) + ...
            delayline(pulse,delay(xdx,idx)*fs,weight(xdx,idx),conf);
    end
end
