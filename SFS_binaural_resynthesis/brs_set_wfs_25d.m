function brs = brs_set_wfs_25d(X,Y,phi,xs,ys,L,src,irs,conf)
%BRS_SET_WFS_25D generates a BRS set for use with the SoundScapeRenderer
%   Usage: brs = brs_set_wfs_25d(X,Y,phi,xs,ys,L,src,irs,conf)
%          brs = brs_set_wfs_25d(X,Y,phi,xs,ys,L,src,irs)
%
%   Input parameters:
%       X,Y     - listener position (m)
%       phi     - listener direction [head orientation] (rad)
%       xs,ys   - virtual source position [ys > Y0 => focused source] (m)
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
%   BRS_SET_WFS_25D(X,Y,phi,xs,ys,L,irs,conf) prepares a BRS set for 
%   a virtual source at [xs ys] for a linear WFS array and the given 
%   listener position.
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
nargmin = 7;
nargmax = 8;
error(nargchk(nargmin,nargmax,nargin));

isargscalar(X,Y,phi,xs,ys);
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
for i = 1:length(angles)
    % Compute BRIR for the desired WFS system
    brs(:,(i-1)*2+1:i*2) = brs_wfs_25d(X,Y,angles(i)+phi,xs,ys,L,src,irs,conf);
end
