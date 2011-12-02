function outsig = auralize_ir_file(irfile,content,conf)
%AURALIZE_IR_FILE Auralize a impulse response wav file with audio content
%   Usage: outsig = auralize_ir_file(brsfile,[file,signal,'content'],conf)
%          outsig = auralize_ir_file(brsfile,[file,signal,'content'])
%
%   Input options:
%       irfile          - file containing impulse response (IR)
%       content         - content file or signal vector to be used for auralisation (mono,
%                         if it contains more than one channel, only the
%                         first will be used).
%                         Also predefined content can be used by applying the
%                         one of the following strings:
%                         'speech', 'noise', 'pinknoise', 'cello', 'castanets'.
%                         Then these contents will be used to auralise the IR.
%                         The corresponding content files are specified in
%                         SFS_config.
%
%   AURALIZE_IR_FILE(irfile,content) convolves the first two channels of
%   the given impulse response file with the given content and returns the resulting
%   outsig. If instead of an explicite content file or vector only a string containig
%   'speech', 'noise', 'pinknoise', 'cello' or 'castanets' is given, the corresponding
%   content file is used.
%
%   see also: auralize_ir, brs_wfs_25d
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input parameters and configuration =================
nargmin = 1:
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargfile(irfile);


%% ===== Computation ====================================================
% Read the IR file
[ir,conf.fs] = wavread(irfile);
% Convolve the file
outsig = auralize_ir(ir,content,conf);
