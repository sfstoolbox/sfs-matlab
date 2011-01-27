function short_ir = shorten_ir(ir,fs,nsamples,conf)
%SHORTEN_HRIR Shortens a HRIR
%   Usage: short_ir = shorten_ir(ir,fs,nsamples,conf)
%          short_ir = shorten_ir(ir,fs,nsamples)
%
%   Input parameters:
%       ir          - two channel IR signal
%       fs          - sampling rate of the target HRIR
%       nsamples    - length of the target HRIR
%       conf        - optional struct containing configuration variables
%                     (see SFS_config for default values)
%
%   Output paramteres:
%       short_ir    - two channel IR signal
%
%   SHORTEN_HRIR(ir,fs,nsamples,conf) shortens a given IR by resampling
%   and applying a hanning window. This is useful e.g. for mobile phones.
%
%   see also: SFS_config, read_irs, ir_intpol
%

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input  parameters ==================================

nargmin = 3;
nargmax = 4;
if nargchk(nargmin,nargmax,nargin)
    error(['Wrong number of args. Usage: short_ir = '],...
        ['shorten_ir(ir,fs,nsamples,conf)']);
end

if nargin<nargmax
    useconfig = true;
elseif ~isstruct(conf)
    error('%s: conf has to be a struct.',upper(mfilename));
else
    useconfig = false;
end


%% ===== Configuration ==================================================

% Load default configuration values
if(useconfig)
    conf = SFS_config;
end

ofs = conf.fs;  % original fs
useplot = conf.useplot;


%% ===== Computation ====================================================

% Resample HRIR
resamp_ir(:,1) = resample(ir(:,1),fs,ofs);
resamp_ir(:,2) = resample(ir(:,2),fs,ofs);

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

