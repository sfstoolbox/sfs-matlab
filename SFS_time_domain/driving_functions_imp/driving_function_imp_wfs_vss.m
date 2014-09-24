function d = driving_function_imp_wfs_vss(x0,xv,dv, conf)
%DRIVING_FUNCTION_MONO_WFS_VSS returns the driving signal D for a virtual
%secondary source distribution
%
%   Usage: D = driving_function_mono_wfs_vss(x0,xv,Dv,f,conf)
%
%   Input parameters:
%       x0          - position, direction, and weights of the real secondary
%                     sources / m [nx7]
%       xv          - position, direction, and weights of the virtual secondary
%                     sources / m [mx7]
%       Dv          - driving functions of virtual secondary sources [mx1]
%       f           - frequency of the monochromatic source / Hz
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   References:
%       S. Spors (2010) - "Local Sound Field Synthesis by Virtual Secondary
%                          Sources", 40th AES
%
%   see also: driving_function_mono_wfs, driving_function_mono_wfs_fs

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
vsstype = conf.localsfs.vss.type;
fs = conf.fs;
N = conf.N;

%% ===== Computation ====================================================
% Get driving signals
if strcmp('ps',vsstype)
    % === Focussed Point Sink ===========================================
    conf.driving_functions = 'default';
elseif strcmp('ls',vsstype)
    % === Focussed Line Sink ============================================
    % Driving signal
    conf.driving_functions = 'default';
else
    error('%s: %s is not a known source type.',upper(mfilename), vsstype);
end

% Get Drivings Signal real secondary sources
N0 = size(x0,1);

d = zeros(N, N0);

idx = 0;
for xvi = xv'
  idx = idx + 1;
  [x0s, xdx] = secondary_source_selection(x0,xvi(1:6)','fs');
  if ~isempty(x0s) && xvi(7) > 0
    % wfs pre equalization of virtual secondary sources driving signal
    pulse = repmat(wfs_preequalization(dv(:,idx), conf), 1, size(x0s,1));
    % pulse = repmat(dv(:,idx), 1, size(x0s,1));
    
    [~, delay, weight] = ...
      driving_function_imp_wfs(x0s,xvi(1:6)','fs',conf);

    x0s = secondary_source_tapering(x0s,conf);

    weight = weight.*x0s(:,7).*xvi(7);
    
    % Shift and weight prototype driving function
    d(:, xdx) = d(:, xdx) + delayline(pulse, delay*fs, weight, conf);
  end
end


