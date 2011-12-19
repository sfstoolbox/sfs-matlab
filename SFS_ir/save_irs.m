function save_irs(irs,outfile)
%SAVE_IRS Save a HRIR/BRIR dataset
%   Usage: save_irs(irs,irsfile)
%
%   Input parameters:
%       irs     - irs struct
%       irsfile - filename of irs mat file
%
%   SAVE_IRS(irs,irsfile) saves a IR dataset as a struct containing the format 
%   specific fields. For a description of the mat format for the IR datasets, 
%   see IR_format.txt.
%   
%   see also: read_irs, get_ir, intpol_ir, dummy_irs, new_irs
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
check_irs(irs);
isargchar(outfile);


%% ===== Save IR file ===================================================

% Save as mat file
% If -v7 doesn't worj use -v6, but note this won't use any compression of your
% data
save('-v7',outfile,'irs');
