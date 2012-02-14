function isargnumeric(varargin)
%ISARGNUMERIC tests if the given arg is numeric and returns an error otherwise
%   Usage: isargnumeric(arg1,arg2,...)
%
%   Input options:
%       args        - list of args
%
%   ISARGNUMRIC(args) tests if all given args are numeric and returns
%   an error otherwise.
%
%   see also: isargscalar

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking for numeric value ======================================
for ii = 1:nargin
    if ~isnumeric(varargin{ii})
        error('%s need to be numeric.',inputname(ii));
    end
end
