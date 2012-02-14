function isargsecondarysource(varargin)
%ISARGSECONDARYSOURCE tests if the given arg is a matrix containing secondary
%   source positions and directions and return an error otherwise
%
%   Usage: isargsecondarysource(args)
%
%   Input options:
%       args        - list of args
%
%   ISARGSECONDARYSOURCE(args) tests if all given args are a matrix
%   containing secondary source positions and directions and returns 
%   an error otherwise.
%
%   see also: isargmatrix
%

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking for vector =============================================
for ii = 1:nargin
    if ~isnumeric(varargin{ii}) || ~ismatrix(varargin{ii}) || ...
        size(varargin{ii},2)~=6
        error(['%s need to be a nx6 matrix containing the secondary ', ...
            'sources positions and directions.'],inputname(ii));
    end
end
