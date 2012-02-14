function isargfile(varargin)
%ISARGFILE tests if the given arg is a file and returns an error otherwise
%   Usage: isargfile(arg1,arg2,...)
%
%   Input options:
%       args        - list of args
%
%   ISARGFILE(args) tests if all given args are a file and returns 
%   an error otherwise.
%
%   see also: isargdir

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking for file =============================================
for ii = 1:nargin
    if ~ischar(varargin{ii}) || ~exist(varargin{ii},'file')
        error('%s need to be a file.',inputname(ii));
    end
end
