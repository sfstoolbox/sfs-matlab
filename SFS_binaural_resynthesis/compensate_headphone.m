function brir = compensate_headphone(brir,conf)
%COMPENSATE_HEADPHONE applies a headphone compensation to the BRIR
%   Usage: brir = compensate_headphone(brir,conf)
%          brir = compensate_headphone(brir)
%
%   Input parameters:
%       brir    - BRIR to which the compensation should be applied
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output:
%       brir    - BRIR which is compensated for the given headphone
%
%   COMPENSATE_HEADPHONE(brir,conf) applies a headphone compensation to the
%   given BRIR. Which headphone compensation it should use is mentioned in the
%   conf struct.
%   The compensation is only applied, if the conf.usehcomp value is not false.
%
%   see also: SFS_config, brs_wfs_25d

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(brir);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
usehcomp = conf.usehcomp;
hcomplfile = conf.hcomplfile;
hcomprfile = conf.hcomprfile;


%% ===== Computation =====================================================
if(usehcomp)
    % Read headphone compensation filter
    hcompl = wavread(hcomplfile);
    hcompr = wavread(hcomprfile);
    hcomp = [hcompl hcompr];
    % Apply filter
    brir(:,1) = conv(hcomp(:,1),brir(1:end-length(hcomp)+1,1));
    brir(:,2) = conv(hcomp(:,2),brir(1:end-length(hcomp)+1,2));
end
