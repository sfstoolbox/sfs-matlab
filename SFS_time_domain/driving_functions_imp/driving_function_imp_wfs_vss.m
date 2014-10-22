function d = driving_function_imp_wfs_vss(x0,xv,dv,conf)
%DRIVING_FUNCTION_IMP_WFS_VSS returns the driving signal d for a given set of
%virtual secondary sources and the corresponding dricing signals
%
%   Usage: d = driving_function_imp_wfs_vss(x0,xv,dv,conf)
%
%   Input parameters:
%       x0          - position, direction, and weights of the real secondary
%                     sources / m [nx7]
%       xv          - position, direction, and weights of the virtual secondary
%                     sources / m [mx7]
%       dv          - driving signals of virtual secondary sources [sxm]
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       d           - driving function signal [Sxn]
%
%   References:
%       S. Spors (2010) - "Local Sound Field Synthesis by Virtual Secondary
%                          Sources", 40th AES
%
%   see also: driving_function_imp_localwfs, driving_function_mono_wfs_vss

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************

%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargmatrix(dv);
isargsecondarysource(x0,xv);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end

%% ===== Configuration ==================================================
dimension = conf.dimension;
fs = conf.fs;
N = conf.N;

%% ===== Computation ====================================================
% Apply wfs preequalization filter on each driving signal of the vss'
dv = wfs_preequalization(dv, conf);

% initialize 
Nv = size(xv,1);
N0 = size(x0,1);

d = zeros(N, N0);
delay = inf(N0,Nv);
weight = zeros(N0,Nv);

idx = 1;
for xvi = xv'
  % select active source for one focused source
  [x0s, xdx] = secondary_source_selection(x0,xvi(1:6)','fs');
  if ~isempty(x0s) && xvi(7) > 0
    % focused source position
    xs = repmat(xvi(1:3)',[size(x0s,1) 1]);    
    % delay and weights for single focused source
    [delay(xdx,idx),weight(xdx,idx)] = driving_function_imp_wfs_fs(x0s(:,1:3),x0s(:,4:6),xs,conf);
    
    % optional tapering
    x0s = secondary_source_tapering(x0s,conf);
    % apply secondary sources' tapering and possibly virtual secondary
    % sources' tapering to weighting matrix
    weight(xdx,idx) = weight(xdx,idx).*x0s(:,7).*xvi(7);
  end
  idx = idx + 1;  
end

% Remove delay offset, in order to begin always at t=0 with the first wave front
% at any secondary source
delay = delay - min(delay(:));

% compose impulse responses
for idx=1:Nv
  xdx = weight(:,idx) ~= 0;
  if sum(xdx) > 0    
    % Shift and weight prototype driving function
    pulse = repmat(dv(:,idx), 1, sum(xdx));
    d(:, xdx) = d(:, xdx) + delayline(pulse, delay(xdx,idx)*fs, weight(xdx,idx), conf);
  end
end
