function isargvector(varargin)
%ISARGVECTOR tests if the given arg is a vector and returns an error otherwise
%   Usage: isargvector(args)
%
%   Input options:
%       args        - list of args
%
%   ISARGVECTOR(args) tests if all given args are a vector and returns 
%   an error otherwise.
%
%   see also: isargscalar, isargmatrix

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking for vector =============================================
for ii = 1:nargin
    if ~isvector(varargin{ii})
        error('%s need to be a vector.',inputname(ii));
    end
end
