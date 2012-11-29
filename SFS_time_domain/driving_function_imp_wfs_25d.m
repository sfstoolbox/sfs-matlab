function [weight,delay] = driving_function_imp_wfs_25d(x0,xs,src,conf)
%DRIVING_FUNCTION_IMP_WFS_25D calculates the WFS 2.5D weighting and delaying
%
%   Usage: [weight,delay] = driving_function_imp_wfs_25d(x0,xs,src,[conf]);
%
%   Input parameters:
%       x0      - position  and direction of secondary source (m)
%       xs      - position of virtual source or diirection of plane wave (m)
%       src     - source type of the virtual source
%                     'pw' - plane wave (xs, ys are the direction of the
%                            plane wave in this case)
%                     'ps' - point source
%                     'fs' - focused source
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       weight  - weight (amplitude) of the driving function
%       delay   - delay of the driving function (s)
%
%   DRIVING_FUNCTION_IMP_WFS_25D(x0,xs,src,conf) returns the
%   weighting and delay parameters of the WFS 2.5D driving function for the given
%   source type and position and loudspeaker positions.
%
%   see also: wave_field_imp_wfs_25d, driving_function_mono_wfs_25d

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargsecondarysource(x0)
isargposition(xs);
xs = position_vector(xs);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
% Speed of sound
c = conf.c;
xref = position_vector(conf.xref);


%% ===== Computation =====================================================

% Direction and position of secondary sources
nx0 = x0(4:6);
x0 = x0(1:3);

% Constant amplitude factor
g0 = sqrt(2*pi*norm(xref-x0));

if strcmp('pw',src)
    % === Plane wave ===
    % Direction of plane wave
    nxs = xs / norm(xs);
    % Delay and amplitude weight
    % NOTE: <n_pw,n(x0)> is the same as the cosinus between their angle
    delay = 1/c * nxs*x0';
    weight = 2*g0 * nxs*nx0';

elseif strcmp('ps',src)
    % === Point source ===
    % Distance between loudspeaker and virtual source
    r = norm(x0-xs);
    % Delay and amplitude weight
    delay = r/c;
    weight = g0/(2*pi)*(x0-xs)*nx0'*r^(-3/2);

elseif strcmp('fs',src)
    % === Focused source ===
    % Distance between loudspeaker and virtual source
    r = norm(x0-xs);
    % Delay and amplitude weight
    delay =  -r/c;
    weight = g0/(2*pi)*(x0-xs)*nx0'*r^(-3/2);
else
    error('%s: %s is not a known source type.',upper(mfilename),src);
end
