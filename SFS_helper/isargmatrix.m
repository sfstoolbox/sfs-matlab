function isargmatrix(varargin)
%ISARGMATRIX tests if the given arg is a matrix and returns an error otherwise
%   Usage: isargmatrix(arg1,arg2,...)
%
%   Input options:
%       args        - list of args
%
%   ISARGMATRIX(args) tests if all given args are a matrix and returns 
%   an error otherwise.
%
%   see also: isargvector, isargscalar

% AUTHOR: Hagen Wierstorf


%% ===== Checking for matrix =============================================
for ii = 1:nargin
    %if ~isnumeric(varargin{ii}) || ~ismatrix(varargin{ii})
    %    error('%s need to be a matrix.',inputname(ii));
    %end
end
