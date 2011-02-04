function outsig = auralize_brs_file(brsfile,contentfile)
%AURALIZE_BRS_FILE Auralize a BRS with an audio file
%   Usage: outsig = auralize_brs_file(brsfile)
%          outsig = auralize_brs_file(brsfile,'speech')
%          outsig = auralize_brs_file(brsfile,'cello')
%          outsig = auralize_brs_file(brsfile,'castanets')
%          outsig = auralize_brs_file(brsfile,contentfile)
%
%   Input options:
%       brsfile         - BRS file
%       contentfile     - content file to be used for auralisation (mono).
%                         Also the strings 'speech', 'cello' and
%                         'castanets' are possible, then these contents
%                         will be used to auralise the BRS.
%
%   AURALIZE_BRS_FILE(brsfile,contentfile) convolves the first two channels of
%   the given brsfile with the given contentfile and returns the resulting
%   outsig as a auralisation of the BRS. If instead of an explicite contentfile 
%   only a string containig 'speech', 'cello' or 'castanets' is given, the 
%   corresponding contentfile is used. If no contentfile is given the
%   castanets contentfile is used.
%
%   see also: wfs_brs, auralize_brs
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input parameters  and configuration ================
nargmin = 1:
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargfile({brsfile},{'brsfile'});


%% ===== Configuration ==================================================

% Load default configuration values
conf = SFS_config;

% Auralisation files
speechfile = conf.speechfile;
cellofile = conf.cellofile;
castanetsfile = conf.castanetsfile;
noisefile = conf.noisefile;
pinknoisefile = conf.pinknoisefile;

if nargin<2 || strcmp(contentfile,'castanets')
    contentfile = castanetsfile;
elseif strcmp(contentfile,'speech')
    contentfile = speechfile;
elseif strcmp(contentfile,'cello')
    contentfile = cellofile;
elseif strcmp(contentfile,'noise')
    contentfile = noisefile;
elseif strcmp(contentfile,'pinknoise')
    contentfile = pinknoisefile;
end

if ~exist(contentfile,'file')
    error('%s: contentfile was not found.',upper(mfilename));
end


%% ===== Computation ====================================================

% Read the content file
[content,contentfs,contentnbit] = wavread(contentfile);

% Check if the content vector has the right format
if ~isvector(content)
    error('%s: content has to be a one channel signal.',upper(mfilename));
end

% Read the BRS file
[brs,brsfs,brsnbit] = wavread(brsfile);

% Check if we have to resample
if contentfs ~= brsfs
    content = resample(content,brsfs,contentfs);
end

% Convolve the two
outsig(:,1) = conv(brs(:,1),content);
outsig(:,2) = conv(brs(:,2),content);

% Scale output
outsig = 0.95*outsig/max(abs(outsig(:)));
