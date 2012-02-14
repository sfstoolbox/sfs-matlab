function varargout = row_vector(varargin)
%ROW_VECTOR make a row vector from the given vector
%
%   Usage: [x1,x2,...] = row_vector(x1,x2,...)
%
%   Input parameters:
%       x1,x2,...  - one or more vectors
%
%   Output parameters:
%       x1,x2,...  - input vectors as row vectors
%
%   ROW_VECTOR(x1,x2,...) returns the given vectors x1,x2,... as row
%   vectors.
%
% see also: column_vector
%

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ===================================
isargvector(varargin{:});


%% ===== Computation =====================================================
for ii = 1:nargin
    if size(varargin{ii},2)<size(varargin{ii},1)
        varargout{ii} = varargin{ii}';
    else
        varargout{ii} = varargin{ii};
    end
end
