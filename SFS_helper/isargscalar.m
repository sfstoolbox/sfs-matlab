function isargscalar(varargin)
%ISARGSCALAR tests if the given arg is a scalar and returns an error otherwise
%   Usage: isargscalar(arg1,arg2,...)
%
%   Input options:
%       args        - list of args
%
%   ISARGSCALAR(args) tests if all given args are a scalar and returns 
%   an error otherwise.
%
%   see also: isargpositivescalar

% AUTHOR: Hagen Wierstorf


%% ===== Checking for scalar =============================================
for ii = 1:nargin
    if ~isnumeric(varargin{ii}) || ~isscalar(varargin{ii})
        error('%s need to be a scalar.',inputname(ii));
    end
end
