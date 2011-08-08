function isargcoord(varargin)
%ISARGCOORD tests if the given arg is a point in an euclidian coordinate system
%
%   Usage: isargcoord(args)
%
%   Input options:
%       args        - list of args
%
%   ISARGCOORD(args) tests if all given args are a point in an euclidian
%   coordinate system. This means they have to be a vector with a length of 1,2
%   or 3. Returns an error otherwise.
%
%   see also: isargvector

% AUTHOR: Hagen Wierstorf


%% ===== Checking for vector =============================================
for ii = 1:nargin
    if ~isnumeric(varargin{ii}) || ~isvector(varargin{ii}) || ...
        length(varargin{ii})>3
        error('%s need to be a vector with length <4.',inputname(ii));
    end
end
