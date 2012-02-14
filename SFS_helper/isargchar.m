function isargchar(varargin)
%ISARGCHAR tests if the given arg is a char and returns an error otherwise
%   Usage: isargstruct(arg1,arg2,...)
%
%   Input options:
%       args        - list of args
%
%   ISARGCHAR(args) tests if all given args are a char and returns 
%   an error otherwise.
%
%   see also: isargstruct

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking for struct =============================================
for ii = 1:nargin
    if ~ischar(varargin{ii})
        error('%s need to be a string.',inputname(ii));
    end
end
