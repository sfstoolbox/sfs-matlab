function addirspath(varargin)
%ADDIRSPATH adds directories containing irs files to the path
%   Usage: addirspath(varargin)
%
%   Input parameters:
%       varargin - path or pathe containing irs data sets. 
%                  Default: '~/svn/ir_databases' and '~/svn/measurements'
%
%   ADDIRSPATH(varargin) adds the given directorysand its subdirectories to
%   the path. If no directory is given, '~/svn/ir_databases' and 
%   '~/svn/measurements' are added.

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


% ===== Checking of input parameters ====================================
nargmin = 0;
nargmax = inf;
error(nargchk(nargmin,nargmax,nargin));
if nargin==nargmin
    dirs{1} = '~/svn/ir_databases';
    dirs{2} = '~/svn/measurements';
else
    isargdir(char(varargin));
    dirs = varargin;
end


%% ===== Adding pathes ==================================================
% FIXME: this adds also all .svn subdirectories!
for ii = 1:length(dirs)
    path = genpath(dirs{ii});
    addpath(path,'-end');
end
