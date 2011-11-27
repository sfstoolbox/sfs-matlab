function outsig = auralize_brs(brs,contentfile,conf)
%AURALIZE_BRS Auralizes a BRS with an audio file
%   Usage: outsig = auralize_brs(brs,contentfile,conf)
%          outsig = auralize_brs(brs,contentfile)
%          outsig = auralize_brs(brs,'speech')
%          outsig = auralize_brs(brs,'cello')
%          outsig = auralize_brs(brs,'castanets')
%
%   Input parameters:
%       brs             - input BRS
%       contentfile     - content file to be used for auralisation (mono,
%                         if it contains more than one channel, only the
%                         first will be used).
%                         Also the strings 'speech', 'cello' and
%                         'castanets' are possible, then these contents
%                         will be used to auralise the BRS.
%
%   AURALIZE_BRS(brs,contentfile) convolves the first two channels of the given
%   BRS with the given contentfile and returns the resulting outsig as a
%   auralisation of the BRS. If instead of an explicite contentfile only a
%   string containig 'speech', 'cello' or 'castanets' is given, the
%   corresponding contentfile is used. If no contentfile is given the
%   castanets contentfile is used.
%
%   see also: auralize_brs_file, wfs_brs
%
% FIXME: change this file to auralize_ir and allow also direct signal as
% input

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input parameters and configuration =================
nargmin = 2;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(brs);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================

% Sampling rate
fs = conf.fs;

% Auralisation files
speechfile = conf.speechfile;
cellofile = conf.cellofile;
castanetsfile = conf.castanetsfile;
noisefile = conf.noisefile;
pinknoisefile = conf.pinknoisefile;

if nargin<2 || strcmp(contentfile,'castanets')
    contentfile = castanetsfile;
elseif isvector(contentfile)
    content = contentfile;
    contentfs = conf.fs;
elseif strcmp(contentfile,'speech')
    contentfile = speechfile;
elseif strcmp(contentfile,'cello')
    contentfile = cellofile;
elseif strcmp(contentfile,'noise')
    contentfile = noisefile;
elseif strcmp(contentfile,'pinknoise')
    contentfile = pinknoisefile;
end

%if ~exist(contentfile,'file')
%    error('%s: contentfile was not found.',upper(mfilename));
%end


%% ===== Computation ====================================================

% Read the content file
if ~exist('content','var')
    [content,contentfs,contentnbit] = wavread(contentfile);
end

% Check if the content vector has the right format
if ~isvector(content)
    error('%s: content has to be a one channel signal.',upper(mfilename));
end

% Check if we have to resample
if contentfs~=fs
    content = resample(content,fs,contentfs);
end

% Convolve the two
for ii = 1:size(brs,2)
    outsig(:,ii) = conv(brs(:,ii),content);
end

% Scale output
outsig = 0.95*outsig/max(abs(outsig(:)));
