function addirspath(basepath)
%ADDIRSPATH adds the directory containing irs files to the path
%   Usage: addirspath(basepath)
%
%   Input options:
%       basepath - path where the databases are located. Default:
%                  ~/svn/ir_databases
%
%   ADDIRSPATH(basepath) add basepath and its subdirectories to the path. If
%   basepath is omitted, ~/svn/ir_databases is used.

% AUTHOR: Hagen Wierstorf


% ===== Checking of input parameters ====================================
nargmin = 0;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
if nargin==nargmin
    basepath = '~/svn/ir_databases';
else
    isargdir(basepath);
end


%% ===== Adding pathes ==================================================
% FIXME: this adds also all .svn subdirectories!
path = genpath(basepath);
addpath(path,'-end');
