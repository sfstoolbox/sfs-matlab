function short_ir = shorten_ir(ir,fs,nsamples,conf)
%SHORTEN_IR shortens a IR
%   Usage: short_ir = shorten_ir(ir,fs,nsamples,conf)
%          short_ir = shorten_ir(ir,fs,nsamples)
%
%   Input parameters:
%       ir          - two channel IR signal
%       fs          - sampling rate of the target IR
%       nsamples    - length of the target IR
%       conf        - optional struct containing configuration variables
%                     (see SFS_config for default values)
%
%   Output paramteres:
%       short_ir    - two channel IR signal
%
%   SHORTEN_HRIR(ir,fs,nsamples,conf) shortens a given IR by resampling
%   and applying a hanning window. This is useful e.g. for mobile phones.
%
%   see also: SFS_config, read_irs, intpol_ir
%

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 4;
error(nargchk(nargmin,nargmax,nargin));

isargpositivescalar(fs,nsamples);
if ~isnumeric(ir) || size(ir,2)~=2
    error('%s: ir has to be an IR with samples x 2 size.',upper(mfilename));
end

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================

ofs = conf.fs;  % original fs
useplot = conf.useplot;


%% ===== Computation ====================================================

% Resample HRIR
if ofs~=fs
    resamp_ir(:,1) = resample(ir(:,1),fs,ofs);
    resamp_ir(:,2) = resample(ir(:,2),fs,ofs);
else
    resamp_ir = ir;
end

% Window HRIR
win = hanningwin(ceil(0.15*nsamples),ceil(0.10*nsamples),nsamples).^2;

% Find maximum of resampled HRIR
% Find maximum in each channel and calculate the mean of the index
[a,idx1] = max(abs(resamp_ir(:,1)));
[a,idx2] = max(abs(resamp_ir(:,2)));
idx = round((idx1+idx2)/2);

% Cut the HRIR around the maximum
% Leading zeros before idx
offset = 24;
short_ir(1:nsamples,1) = ...
    resamp_ir(idx-offset:idx+nsamples-offset-1,1) .* win;
short_ir(1:nsamples,2) = ...
    resamp_ir(idx-offset:idx+nsamples-offset-1,2) .* win;


%% ===== Plotting =======================================================

if(useplot)
    figure
    plot(resamp_ir(:,1),'-b'); hold on;
    plot(resamp_ir(:,2),'-r');
    figure
    plot(short_ir(:,1),'-b'); hold on;
    plot(short_ir(:,2),'r-');
end

