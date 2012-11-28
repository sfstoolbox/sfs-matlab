function brs = brs_wfs_25d(X,phi,xs,src,L,irs,conf)
%BRS_WFS_25D generates a BRS set for use with the SoundScapeRenderer
%
%   Usage: brs = brs_wfs_25d(X,phi,xs,src,L,irs,[conf])
%
%   Input parameters:
%       X       - listener position (m)
%       phi     - listener direction [head orientation] (rad)
%       xs      - virtual source position [ys > Y0 => focused source] (m)
%       src     - source type: 'pw' - plane wave
%                              'ps' - point source
%                              'fs' - focused source
%       L       - Length of linear loudspeaker array (m)
%       irs     - IR data set for the second sources
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       brs     - conf.N x 2*nangles matrix containing all brs (2
%                 channels) for every angles of the BRS set
%
%   BRS_WFS_25D(X,phi,xs,src,L,irs,conf) prepares a BRS set for
%   a virtual source at xs for a linear WFS array and the given
%   listener position.
%   One way to use this BRS set is using the SoundScapeRenderer (SSR), see
%   http://www.tu-berlin.de/?id=ssr
%
%   Geometry (for a linear array):
%
%                               y-axis
%                                 ^
%                                 |
%                                 |
%                                 |    (Listener)
%                                 |        O X, phi=-pi/2
%                                 |        |
%                                 |
%                  o xs           |
%             (Virtual Source)    |
%                                 |
%     -------v--v--v--v--v--v--v--v--v--v--v--v--v--v--v------> x-axis
%                                 X0 (Array center)
%            |---      Loudspeaker array length     ---|
%
%   see also: SFS_config, brs_wfs_25d, brs_point_source

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
nargmin = 6;
nargmax = 7;
narginchk(nargmin,nargmax);
[X,xs] = position_vector(X,xs);
isargscalar(phi);
isargpositivescalar(L);
check_irs(irs);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
N = conf.N;                     % Target length of BRIR impulse responses
angles = rad(conf.brsangles);   % Angles for the BRIRs


%% ===== Computation =====================================================
% Initial values
brs = zeros(N,2*length(angles));

% Generate a BRS set for all given angles
for ii = 1:length(angles)
    % Compute BRIR for the desired WFS system
    brs(:,(ii-1)*2+1:ii*2) = ...
        ir_wfs_25d(X,angles(ii)+phi,xs,src,L,irs,conf);
end
