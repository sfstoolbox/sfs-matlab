function isargmatrix(args,argnames)
%ISARGMATRIX tests if the given arg is a matrix and returns an error otherwise
%   Usage: isargmatrix(args,argnames)
%
%   Input options:
%       args        - list of args
%       argnames    - list of corresponding argnames (strings)
%
%   ISARGMATRIX(args,argnames) tests if all given args are a matrix and returns 
%   an error otherwise.
%
%   see also: isargvector, isargscalar

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

%% ===== Checking for matrix =============================================
for ii = 1:length(args)
    if ~isnumeric(args{ii}) || ~ismatrix(args{ii})
        error('%s need to be a matrix.',argnames{ii});
    end
end
