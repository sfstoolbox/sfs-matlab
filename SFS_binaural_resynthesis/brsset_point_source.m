function brs = brsset_point_source(X,phi,xs,irs,conf)
%BRSSET_POINT_SOURCE generates a BRS set for use with the SoundScapeRenderer
%   Usage: brs = brsset_point_source(X,phi,xs,irs,conf)
%          brs = brsset_point_source(X,phi,xs,irs,)
%
%   Input parameters:
%       X       - listener position (m)
%       phi     - listener direction [head orientation] (rad)
%       xs      - source position (m)
%       irs     - IR data set for the second sources
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameters:
%       brs     - conf.N x 2*nangles matrix containing all brs (2
%                 channels) for every angles of the BRS set
%
%   BRSSET_POINT_SOURCE(X,phi,xs,irs,conf) prepares a BRS set for
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
%   see also: SFS_config, brs_point_source, brs_wfs_25d
%

% AUTHOR: Sascha Spors, Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 5;
error(nargchk(nargmin,nargmax,nargin));
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
for i = 1:length(angles)
    % Compute BRIR for a reference (single loudspeaker at xs)
    brs(:,(i-1)*2+1:i*2) = brs_point_source(X,angles(i)+phi,xs,irs,conf);
end
warning('on','SFS:irs_intpol');

