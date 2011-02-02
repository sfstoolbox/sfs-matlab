function irs = read_irs(conf)
%READ_IRS Read a HRIR/BRIR dataset
%   Usage: irs = read_irs(conf)
%          irs = read_irs()
%
%   Input parameters:
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output paramteres:
%       irs   - struct containing
%           .left           - left ear HRIR signal
%           .right          - right ear HRIR signal
%           .angle          - azimuth and elevation angles
%           .r0             - measurement distance of IR dataset (m)
%           .tag            - 'HRIR' or 'BRIR'
%           .description    - description of the IR data set
%
%   READ_IRS(conf) loads a IR dataset as a struct containing the above mentioned
%   fields from the mat files stored in conf.irsfile. For a description of
%   the mat format for the IR datasets, see IR_format.txt.
%   
%   see also: SFS_config, wfs_brs, ref_brs, hrir_intpol, create_irs_mat
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================

if nargchk(0,1,nargin)
    error('Wrong number of args. Usage: irs = read_irs(conf)');
end

if nargin<1
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

irsfile = conf.irsfile;


%% ===== Read IR files ================================================

% Load the mat file
load(irsfile);

% Check irs format
check_irs(irs);
