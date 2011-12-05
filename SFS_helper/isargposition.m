function isargposition(varargin)
%ISARGPOSITION tests if the given arg is a vector and returns an error otherwise
%   Usage: isargposition(args)
%
%   Input options:
%       args        - list of args
%
%   ISARGPOSITION(args) tests if all given args are a position vector and
%   returns an error otherwise.
%
%   see also: isargscalar, isargvector

% AUTHOR: Hagen Wierstorf


%% ===== Checking for vector =============================================
for ii = 1:nargin
    if ~isvector(varargin{ii}) || length(varargin{ii})<2 || ...
        length(varargin{ii})>3
        error('%s need to be a position vector.',inputname(ii));
    end
end
