function isargnegativescalar(varargin)
%ISARGNEGATIVESCALAR tests if the given arg is a negative scalar
%   Usage: isargnegativescalar(arg1,arg2,...)
%
%   Input options:
%       args        - list of args
%
%   ISARGNEGATIVESCALAR(args) tests if all given args are a negative
%   scalar and returns an error otherwise.
%
%   see also: isargscalar, isargpositivescalar

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking for scalar =============================================
for ii = 1:nargin
    if ~isnumeric(varargin{ii}) || ~isscalar(varargin{ii}) || varargin{ii}>0
        error('%s need to be a negative scalar.',inputname(ii));
    end
end
