function D = driving_function_mono_wfs_vss(x0,xv,Dv,f,conf)
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
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
isargvector(Dv);
isargpositivescalar(f);
isargsecondarysource(x0,xv);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end

%% ===== Configuration ==================================================
vsstype = conf.localsfs.vss.type;

%% ===== Computation ====================================================
% Get driving signals
if strcmp('ps',vsstype)
    % === Focussed Point Sink ===========================================
    conf.driving_functions = 'default';
elseif strcmp('ls',vsstype)
    % === Focussed Line Sink ============================================
    % Driving signal
    conf.driving_functions = 'line_sink';
else
    error('%s: %s is not a known source type.',upper(mfilename), vsstype);
end

% Get Drivings Signal real secondary sources
Ns = size(xv,1);
N0 = size(x0,1);

Dmatrix = zeros(N0,Ns);

for idx=1:Ns
  [xtmp, xdx] = secondary_source_selection(x0,xv(idx,1:6),'fs');
  if (~isempty(xtmp))
    xtmp = secondary_source_tapering(xtmp,conf);  
    Dmatrix(xdx,idx) = driving_function_mono_wfs(xtmp,xv(idx,1:3),'fs',f,conf);
    Dmatrix(xdx,idx) = Dmatrix(xdx,idx).*xtmp(:,7);
  end
end

D = Dmatrix*(Dv.*xv(:,7));