function irs = create_android_irs_mat(irset,nsamples,fs,conf)
%CREATE_ANDROID_IRS_MAT generates an IR mat file for the Android phone
%   Usage: short_irs = create_android_irs_mat(irset)
%
%   Input parameters:
%       irset       - Name of the IR set to be shortend. The name is used
%                     in this file to identify the code to run and is used
%                     as a name for the output dir containing the files for
%                     the Android.
%       nsamples    - length of the desired IR set in samples
%       fs          - sampling rate of the desired IR set
%
%   Output paramteres:
%       irs   - struct containing
%           .left           - left ear IR signal
%           .right          - right ear IR signal
%           .angle          - [azimuth; elevation] angles (rad)
%           .r0             - measurement distance of IR dataset (m)
%           .tag            - 'HRIR' or 'BRIR'
%           .description    - description of the IR data set
%
%   CREATE_ANDROID_IRS_MAT(irset) reads a IR dataset from wav file etc.,
%   shortens it by downsampling and clipping, and stores it in a mat file
%   in the format needed by SFS. If you have your own IR dataset you will
%   add a section to this file in order to read it and store it in the
%   desired format.
%   This function is mainly to check the shortening algorithm for the
%   Andoird IR data sets.
%
%   IR sets, that are supported at the moment are (given by irset name):
%
%       'FABIAN_pinta_anechoic'
%           HRIR data set measured by Alexander Lindau in Pinta using the FABIAN
%           manikin. The set contains only data for the horizontal plane with an
%           azimuth resolution of 1.25°
%
%       'Wittek_Studio_echoic'
%           BRIR data set measured by Robert Wittek in Studio (Bochum) with
%           the ... manikin. The set contains only data for the horizontal
%           plane with an azimuth resolution of 1° (NOTE: if you need a
%           finer resolution have a look at the original data, there Wittek
%           uses a resolution of just 0.1°)
%
%
%   see also: create_irs_mat, create_android_irs
%

% AUTHOR: Hagen Wierstorf, Sascha Spors


%% ===== Checking of input  parameters ==================================

if nargchk(3,4,nargin)
    error('Wrong number of args. Usage: short_irs = create_android_irs_mat(irset,nsamples,fs)');
end
if ~ischar(irset)
    error('%s: irset has to be a string!',upper(mfilename));
end
if ~isnumeric(nsamples) || ~isscalar(nsamples) || nsamples<1
    error('%s: nsamples has to be a positive scalar.',upper(mfilename));
end
if ~isnumeric(fs) || ~isscalar(fs) || fs<1
    error('%s: fs has to be a positive scalar.',upper(mfilename));
end
if nargin<4
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


%% ===== Variables ======================================================

% Target angles
angles = (-180:1:180)/180*pi;


%% ===== Computation ====================================================

% Choose the desired IR set
if strcmp(irset,'FABIAN_pinta_anechoic')
    % Read IRs
    conf.irsfile = 'D:/data/measurements/HRIRs/FABIAN_pinta_anechoic.mat';
    outdir = 'D:/data/measurements/HRIRs';
    irs = read_irs(conf);
elseif strcmp(irset,'Wittek_Studio_echoic')
    % Read IRs
    conf.irsfile = 'D:/data/measurements/BRIRs/Wittek_Studio_echoic.mat';
    outdir = 'D:/data/measurements/BRIRs';
    irs = read_irs(conf);
end

% Generate a set of short IRs
short_irs = irs;
short_irs.description = [irs.description,...
    ' Shortend for the Android phone.'];
short_irs.left = zeros(nsamples,length(angles));
short_irs.right = zeros(nsamples,length(angles));
for ii = 1:length(angles)
    ir = get_ir(irs,angles(ii));
    tmp = shorten_ir(ir,fs,nsamples,conf);
    short_irs.left(:,ii) = tmp(:,1);
    short_irs.right(:,ii) = tmp(:,2);
    short_irs.angle(:,ii) = [angles(ii); 0];    
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