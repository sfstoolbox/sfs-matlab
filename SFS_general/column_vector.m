function varargout = column_vector(varargin)
%COLUMN_VECTOR make a column vector from the given vector
%   Usage: [x1,x2,...] = column_vector(x1,x2,...)
%
%   Input parameters:
%       x1,x2,...  - one or more vectors
%
%   Output parameters:
%       x1,x2,...  - input vectors as column vectors
%
%   COLUMN_VECTOR(x1,x2,...) returns the given vectors x1,x2,... as column
%   vectors.
%
% see also: row_vector
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ===================================
isargvector(varargin{:});


%% ===== Computation =====================================================
for ii = 1:nargin
    if size(varargin{ii},2)>size(varargin{ii},1)
        varargout{ii} = varargin{ii}';
    else
        varargout{ii} = varargin{ii};
    end
end
