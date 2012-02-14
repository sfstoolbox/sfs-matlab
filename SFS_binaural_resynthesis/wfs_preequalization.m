function ir = wfs_preequalization(ir,conf)
%WFS_PREEQUALIZATION applies a pre-equalization filter for WFS
%   Usage: ir = wfs_prefilter(ir,conf)
%          ir = wfs_prefilter(ir)

%
%   Input parameters:
%       ir      - IR to which the pre-equalization filter should be applied
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output:
%       ir      - IR with applied pre-equalization 
%
%   WFS_PREEQUALIZATION(ir,conf) applies the pre-equalization filter for
%   Wave Field Synthesis to the given impulse response.
%
%   see also: wfs_prefilter, SFS_config, brs_wfs_25d

% AUTHOR: Sascha Spors, Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(ir);
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
for ii = 1:size(ir,2)
    ir(:,ii) = conv(hpre,ir(1:end-length(hpre)+1,ii));
end
