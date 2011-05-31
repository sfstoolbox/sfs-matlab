function brir = wfs_preequalization(brir,conf)
%WFS_PREEQUALIZATION applies a pre-equalization filter for WFS
%   Usage: brir = wfs_prefilter(brir,conf)
%          brir = wfs_prefilter(brir)

%
%   Input parameters:
%       brir    - BRIR to which the pre-equalization filter should be applied
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output:
%       brir    - BRIR with applied pre-equalization 
%
%   WFS_PREEQUALIZATION(brir,conf) applies the pre-equalization filter for
%   Wave Field Synthesis to the given BRIR.
%
%   see also: wfs_prefilter, SFS_config, brs_wfs_25d

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(brir);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
usehpre = conf.usehpre;


%% ===== Computation =====================================================
% Check if we should procide
if ~usehpre
    return;
end
% Get the filter
hpre = wfs_prefilter(conf);
% Apply the filter
brir(:,1) = conv(hpre,brir(1:end-length(hpre)+1,1));
brir(:,2) = conv(hpre,brir(1:end-length(hpre)+1,2));
