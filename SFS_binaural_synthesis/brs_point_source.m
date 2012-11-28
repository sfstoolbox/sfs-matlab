function brs = brs_point_source(X,phi,xs,irs,conf)
%BRS_POINT_SOURCE generates a BRS set for use with the SoundScapeRenderer
%
%   Usage: brs = brs_point_source(X,phi,xs,irs,[conf])
%
%   Input parameters:
%       X       - listener position (m)
%       phi     - listener direction [head orientation] (rad)
%       xs      - source position (m)
%       irs     - IR data set for the second sources
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       brs     - conf.N x 2*nangles matrix containing all brs (2
%                 channels) for every angles of the BRS set
%
%   BRS_POINT_SOURCE(X,phi,xs,irs,conf) prepares a BRS set for
%   a reference source (single point source) for the given listener
%   position.
%   One way to use this BRS set is using the SoundScapeRenderer (SSR), see
%   http://www.tu-berlin.de/?id=ssr
%
%   Geometry:
%
%                                 y-axis
%                                   ^
%                                   |
%                                   |
%                                   |
%                                   |    listener
%                                   |       O X, phi=-pi/2
%                                   |       |
%               source              |
%                 o xs              |
%                                   |
%                                   |
%       ----------------------------|---------------------------> x-axis
%
%   see also: SFS_config, ir_point_source, brs_wfs_25d

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
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
[X,xs] = position_vector(X,xs);
isargscalar(phi);
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
warning('off','SFS:irs_intpol');
for ii=1:length(angles)
    % Compute IR for a reference (single loudspeaker at xs)
    brs(:,(ii-1)*2+1:ii*2) = ir_point_source(X,angles(ii)+phi,xs,irs,conf);
end
warning('on','SFS:irs_intpol');
