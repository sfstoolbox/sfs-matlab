function irs = read_irs(irsfile)
%READ_IRS Read a HRIR/BRIR dataset
%   Usage: irs = read_irs(irsfile)
%
%   Input parameters:
%       irsfile - filename of irs mat file
%
%   Output paramteres:
%       irs   - irs struct. For details on the containing fields have a look at
%               the IR_format.txt file.
%
%   READ_IRS(irsfile) loads a IR dataset as a struct containing the format 
%   specific fields. For a description of the mat format for the IR datasets, 
%   see IR_format.txt.
%   
%   see also: get_ir, intpol_ir, dummy_irs, new_irs, wfs_brs, ref_brs
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
if ~ischar(irsfile) || ~exist(irsfile,'file')
    error('%s: irsfile has to be an existing and valid file.',upper(mfilename));
end


%% ===== Read IR files ================================================

% Load the mat file
load(irsfile);

% Check irs format
check_irs(irs);
