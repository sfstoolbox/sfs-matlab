function outsig = auralize_ir(ir,content,usenorm,conf)
%AURALIZE_IR auralizes an impulse response with an audio file/signal
%
%   Usage: outsig = auralize_ir(ir,[file,sig,'content'],normalize,conf)
%          outsig = auralize_ir(ir,[file,sig,'content'],normalize)
%          outsig = auralize_ir(ir,[file,sig,'content'])
%
%   Input parameters:
%       ir        - impulse response (IR). Also an binaural room scanning
%                   (BRS) matrix can be auralized, then only the two first
%                   channels will be used.
%       content   - content file or signal vector to be used for auralisation
%                   (mono, if it contains more than one channel, only the
%                   first will be used).
%                   Also predefined content can be used by applying the
%                   one of the following strings:
%                   'speech', 'noise', 'pinknoise', 'cello', 'castanets'.
%                   Then these contents will be used to auralise the IR.
%                   The corresponding content files are specified in
%                   SFS_config.
%       normalize - normalize the signal (1 or 0), default: 1
%       conf      - optional struct containing configuration variables (see
%                   SFS_config for default values)
%
%   AURALIZE_IR(ir,content,normalize) convolves the first two channels of the
%   given IR with the given content and returns the resulting outsig. If
%   instead of an explicite content file or vector only a string containig
%   'speech', 'noise', 'pinknoise', 'cello' or 'castanets' is given, the
%   corresponding content file is used.
%
%   see also: auralize_ir_file, brs_wfs_25d, brs_point_source

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input parameters and configuration =================
nargmin = 2;
nargmax = 4;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(ir);

if nargin<3
    usenorm = 1;
    conf = SFS_config;
elseif nargin<nargmax
    conf = SFS_config;
end
isargscalar(usenorm);
isargstruct(conf);


%% ===== Configuration ==================================================
% Sampling rate
fs = conf.fs;
% Auralisation files
speechfile = conf.speechfile;
cellofile = conf.cellofile;
castanetsfile = conf.castanetsfile;
noisefile = conf.noisefile;
pinknoisefile = conf.pinknoisefile;


%% ===== Get the right content ==========================================
if isvector(content) & isnumeric(content)
    contentfs = conf.fs;
else
    if strcmp(content,'castanets')
        contentfile = castanetsfile;
    elseif strcmp(content,'speech')
        contentfile = speechfile;
    elseif strcmp(content,'cello')
        contentfile = cellofile;
    elseif strcmp(content,'noise')
        contentfile = noisefile;
    elseif strcmp(content,'pinknoise')
        contentfile = pinknoisefile;
    elseif ~exist(content,'file')
        error('%s: %s file was not found.',upper(mfilename),content);
    else
        contentfile = content;
    end
    % Read the content file
    [content,contentfs] = wavread(contentfile);
end


%% ===== Convolution of the IR with content =============================
% Check if we have to resample
if contentfs~=fs
    content = resample(content,fs,contentfs);
end
% Check if the content is only one channel and fix it otherwise
if min(size(content))~=1
    if size(content,1)<size(content,2)
        content = content(1,:);
    else
        content = content(:,1);
    end
end
% Convolve the two
for ii = 1:2
    outsig(:,ii) = conv(ir(:,ii),content);
end
% Scale output
if(usenorm)
    outsig = norm_signal(outsig);
end
