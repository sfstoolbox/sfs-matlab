function isargnegativescalar(args,argnames)
%ISARGNEGATIVESCALAR tests if the given arg is a negative scalar
%   Usage: isargnegativescalar(args,argnames)
%
%   Input options:
%       args        - list of args
%       argnames    - list of corresponding argnames (strings)
%
%   ISARGNEGATIVESCALAR(args,argnames) tests if all given args are a negative
%   scalar and returns an error otherwise.
%
%   see also: isargscalar, isargpositivescalar

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
if ~iscell(args)
    error('%s: args need to be a cell array.',upper(mfilename));
end
if ~iscell(argnames) || ~iscellstr(argnames)
    error('%s: argnames need to be a cell array containing only string.',...
        upper(mfilename));
end
if length(args)~=length(argnames)
    error('%s: args and argnames must have the same length.',upper(mfilename));
end

%% ===== Checking for scalar =============================================
for ii = 1:length(args)
    if ~isnumeric(args{ii}) || ~isscalar(args{ii}) || args{ii}>0
        error('%s need to be a negative scalar.',argnames{ii});
    end
end
