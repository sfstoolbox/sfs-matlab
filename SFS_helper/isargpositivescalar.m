function isargpositivescalar(varargin)
%ISARGPOSITIVESCALAR tests if the given arg is a positive scalar
%   Usage: isargpositivescalar(arg1,arg2,...)
%
%   Input options:
%       args        - list of args
%
%   ISARGPOSITIVESCALAR(args) tests if all given args are a positive
%   scalar and returns an error otherwise.
%
%   see also: isargscalar, isargnegativescalar

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking for scalar =============================================
for ii = 1:nargin
    if ~isnumeric(varargin{ii}) || ~isscalar(varargin{ii}) || varargin{ii}<0
        error('%s need to be a positive scalar.',inputname(ii));
    end
end
