function brs = brsset_wfs_25d(X,phi,xs,L,src,irs,conf)
%BRS_SET_WFS_25D generates a BRS set for use with the SoundScapeRenderer
%   Usage: brs = brsset_wfs_25d(X,phi,xs,L,src,irs,conf)
%          brs = brsset_wfs_25d(X,phi,xs,L,src,irs)
%
%   Input parameters:
%       X       - listener position (m)
%       phi     - listener direction [head orientation] (rad)
%       xs      - virtual source position [ys > Y0 => focused source] (m)
%       L       - Length of linear loudspeaker array (m)
%       src     - source type: 'pw' - plane wave
%                              'ps' - point source
%                              'fs' - focused source
%       irs     - IR data set for the second sources
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameters:
%       brs     - conf.N x 2*nangles matrix containing all brs (2
%                 channels) for every angles of the BRS set
%
%   BRSSET_WFS_25D(X,phi,xs,L,irs,conf) prepares a BRS set for
%   a virtual source at xs for a linear WFS array and the given
%   listener position.
%   One way to use this BRS set is using the SoundScapeRenderer (SSR), see
%   http://www.tu-berlin.de/?id=ssr
%
%   Geometry:
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

% AUTHOR: Sascha Spors, Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 6;
error(nargchk(nargmin,nargmax,nargin));
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
        brs_wfs_25d(X,angles(ii)+phi,xs,L,src,irs,conf);
end
