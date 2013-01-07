function [weight,delay] = driving_function_imp_wfs_3d(x0,xs,src,conf)
%DRIVING_FUNCTION_IMP_WFS_3D calculates the WFS 3D weighting and delaying
%
%   Usage: [weight,delay] = driving_function_imp_wfs_3d(x0,xs,src,[conf]);
%
%   Input parameters:
%       x0      - position  and direction of secondary source (m)
%       xs      - position of virtual source or diirection of plane wave (m)
%       src     - source type of the virtual source
%                     'pw3D' - plane wave (xs, ys, zs are the direction of the
%                            plane wave in this case)
%                                         
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       weight  - weight (amplitude) of the driving function
%       delay   - delay of the driving function (s)
%
%   DRIVING_FUNCTION_IMP_WFS_3D(x0,xs,src,conf) returns the
%   weighting and delay parameters of the WFS 3D driving function for the given
%   source type and position and loudspeaker positions.
%
%   see also: ....   (wave_field_imp_wfs_25d,
%   driving_function_mono_wfs_25d)

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

% AUTHOR: Hagen Wierstorf
% $LastChangedDate: 2012-04-24 17:54:07 +0200 (Tue, 24 Apr 2012) $
% $LastChangedRevision: 701 $
% $LastChangedBy: wierstorf.hagen $


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 4;
error(nargchk(nargmin,nargmax,nargin));
if strcmp('spherical',conf.array)
    isargsecondarysourceforsphere(x0);
else
    isargsecondarysource(x0)
end
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
% Check also the activity of the used loudspeaker.

% Weights if needed
equallyPointsWeights = x0(7);
surfaceWeights = x0(8);

% Direction and position of active secondary sources
nx0 = x0(4:6);
x0 = x0(1:3);
    
if strcmp('pw',src)
% === Plane wave ===
% Direction of plane wave
nxs = xs / norm(xs);
% Delay and amplitude weight
delay = 1/c * nxs*x0';
weight = 2 .* nxs*nx0'.*equallyPointsWeights.*surfaceWeights.*10; % *10 because the amplitude is to low to see anything in the plot
        
else
        error('%s: %s is not a known source type.',upper(mfilename),src);
end


end