function isargdir(varargin)
%ISARGDIR tests if the given arg is a directory
%   Usage: isargdirectory(arg1,arg2,...)
%
%   Input options:
%       args        - variable number of args
%
%   ISARGDIR(args) tests if all given args are a directory and
%   returns an error otherwise.
%
%   see also: isargfile

% AUTHOR: Hagen Wierstorf


%% ===== Checking for directory ==========================================
for ii = 1:nargin
    if ~ischar(varargin{ii}) || ~exist(varargin{ii},'dir')
        error('%s need to be a directory.',inputname(ii));
    end
end
