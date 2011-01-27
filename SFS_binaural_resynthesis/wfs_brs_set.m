function brs = wfs_brs_set(X,Y,phi,xs,ys,L,irs,conf)
%WFS_BRS_SET generates a BRS set for the SoundScapeRenderer
%   Usage: brs = wfs_brs_set(X,Y,phi,xs,ys,L,irs,outdir,conf)
%          brs = wfs_brs_set(X,Y,phi,xs,ys,L,irs,outdir)
%
%   Input parameters:
%       X,Y     - listener position (m)
%       phi     - listener direction [head orientation] (°)
%       xs,ys   - virtual source position [ys > Y0 => focused source] (m)
%       L       - Length of linear loudspeaker array (m)
%       irs     - IR data set for the second sources
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameters:
%       brs     - conf.N x 2*nangles matrix containing all brs (2
%                 channels) for every angles of the BRS set
%
%   WFS_BRS_SET(X,Y,phi,xs,ys,L,irs,conf) prepares a BRS set for 
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

if nargchk(7,8,nargin)
    error(['Wrong number of args.',...
           'Usage: brs = wfs_brs_set(X,Y,phi,xs,ys,L,irs,conf)']);
end

if ~isnumeric(X) || ~isscalar(X)
    error('%s: X has to be a scalar!',upper(mfilename));
end

if ~isnumeric(Y) || ~isscalar(Y)
    error('%s: Y has to be a scalar!',upper(mfilename));
end

if ~isnumeric(phi) || ~isscalar(phi)
    error('%s: phi has to be a scalar!',upper(mfilename));
end

if ~isnumeric(xs) || ~isscalar(xs)
    error('%s: xs has to be a scalar!',upper(mfilename));
end

if ~isnumeric(ys) || ~isscalar(ys)
    error('%s: ys has to be a scalar!',upper(mfilename));
end

if ~isnumeric(L) || ~isscalar(L) || L<0
    error('%s: L has to be a positive scalar!',upper(mfilename));
end

if ~isstruct(irs)
    error('%s: irs has to be a struct!',upper(mfilename));
end

if nargin<8
    useconfig = true;
elseif ~isstruct(conf)
    error('%s: conf has to be a struct.',upper(mfilename));
else
    useconfig = false;
end


%% ===== Configuration ===================================================

% Load default configuration values
if(useconfig)
    conf = SFS_config;
end

N = conf.N;                 % Target length of BRIR impulse responses
angles = conf.brsangles;    % Angles for the BRIRs


%% ===== Computation =====================================================

% Initial values
brs = zeros(N,2*length(angles));

% Generate a BRIR set for all given angles
for i = 1:length(angles)
    % Compute BRS for a linear WFS system
    brs(:,(i-1)*2+1:i*2) = wfs_brs(X,Y,angles(i)+phi,xs,ys,L,irs,conf);
end

