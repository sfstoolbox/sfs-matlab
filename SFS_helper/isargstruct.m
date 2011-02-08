function isargstruct(varargin)
%ISARGSTRUCT tests if the given arg is a struct and returns an error otherwise
%   Usage: isargstruct(arg1,arg2,...)
%
%   Input options:
%       args        - list of args
%
%   ISARGSTRUCT(args) tests if all given args are a struct and returns 
%   an error otherwise.
%
%   see also: isargchar

% AUTHOR: Hagen Wierstorf


%% ===== Checking for struct =============================================
for ii = 1:nargin
    if ~isstruct(varargin{ii})
        error('%s need to be a struct.',inputname(ii));
    end
end
