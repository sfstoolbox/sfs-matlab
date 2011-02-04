function isargstruct(args,argnames)
%ISARGSTRUCT tests if the given arg is a struct and returns an error otherwise
%   Usage: isargstruct(args,argnames)
%
%   Input options:
%       args        - list of args
%       argnames    - list of corresponding argnames (strings)
%
%   ISARGSTRUCT(args,argnames) tests if all given args are a struct and returns 
%   an error otherwise.
%
%   see also: isargchar

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

%% ===== Checking for struct =============================================
for ii = 1:length(args)
    if ~isstruct(args{ii})
        error('%s need to be a struct.',argnames{ii});
    end
end
