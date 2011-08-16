function brs = brsset_point_source(X,Y,phi,xs,ys,irs,conf)
%BRSSET_POINT_SOURCE generates a BRS set for use with the SoundScapeRenderer
%   Usage: brs = brsset_wfs_25d(X,Y,phi,xs,ys,irs,conf)
%          brs = brsset_wfs_25d(X,Y,phi,xs,ys,irs,)
%
%   Input parameters:
%       X,Y     - listener position (m)
%       phi     - listener direction [head orientation] (rad)
%       xs,ys   - source position (m)
%       irs     - IR data set for the second sources
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameters:
%       brs     - conf.N x 2*nangles matrix containing all brs (2
%                 channels) for every angles of the BRS set
%
%   BRSSET_POINT_SOURCE(X,Y,phi,xs,ys,irs,conf) prepares a BRS set for
%   a reference source (single point source) for the given listener
%   position.
%   One way to use this BRS set is using the SoundScapeRenderer (SSR), see
%   http://www.tu-berlin.de/?id=ssr
%
%   Geometry:
%              |---      Loudspeaker array length     ---|
%    x-axis                      [X0 Y0] (Array center)
%       <------^--^--^--^--^--^--^--^--^--^--^--^--^--^--^-------
%                                   |
%                 x [xs ys]         |
%           (Single Source)         |
%                                   |        |
%                                   |        O [X Y], phi
%                                   |    (Listener)
%                                   |
%                                   |
%                                   |
%                                   |
%                                   v y-axis
%
%   see also: SFS_config, wfs_brs, ref_brs_set

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 6;
nargmax = 7;
error(nargchk(nargmin,nargmax,nargin));

isargscalar(X,Y,phi,xs,ys);
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
for i = 1:length(angles)
    % Compute BRIR for a reference (single loudspeaker at [xs,ys])
    brs(:,(i-1)*2+1:i*2) = brs_point_source(X,Y,angles(i)+phi,xs,ys,irs,conf);
end
