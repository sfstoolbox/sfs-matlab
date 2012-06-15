function irs = create_android_irs_mat(irs,nsamples,fs,conf)
%CREATE_ANDROID_IRS_MAT generates an IR mat file for the Android phone
%   
%   Usage: irs = create_android_irs_mat(irs,nsamples,fs,[conf])
%
%   Input parameters:
%       irs         - irs set to shorten
%       nsamples    - length of the desired IR set in samples
%       fs          - sampling rate of the desired IR set
%
%   Output paramteres:
%       irs         - shortened and resampled irs set
%
%   CREATE_ANDROID_IRS_MAT(irs,nsamples,fs) shortens the given irs set
%   by downsampling and clipping.
%
%   see also: create_irs_mat, check_irs
%

% AUTHOR: Hagen Wierstorf, Sascha Spors
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 4;
error(nargchk(nargmin,nargmax,nargin));
check_irs(irs);
isargpositivescalar(nsamples,fs);

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Computation ====================================================

% Generate a set of short IRs
angles = irs.apparent_azimuth;
short_irs = irs;
short_irs.description = [irs.description,...
    ' Shortend for the Android phone.'];
short_irs.left = zeros(nsamples,length(angles));
short_irs.right = zeros(nsamples,length(angles));
for ii = 1:length(angles)
    ir = get_ir(irs,angles(ii));
    tmp = reduce_ir(ir,fs,nsamples,conf);
    short_irs.left(:,ii) = tmp(:,1);
    short_irs.right(:,ii) = tmp(:,2);
end

irs = short_irs;

% Normalize the short IR set
maxval = max(max(abs([irs.left(:) irs.right(:)])));
irs.left = irs.left ./ maxval;
irs.right = irs.right ./ maxval;


%% ===== Write the HRIR mat file =========================================
if(0)
outfile = sprintf('%s/android_%i_%s.mat',outdir,nsamples,irset);
save(outfile,'irs');
end
