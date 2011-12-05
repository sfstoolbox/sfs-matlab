function varargout = position_vector(varargin)
%POSITION_VECTOR returns 3x1 vectors
%
%   Usage: [x1,x2,...] = position_vector(x1,x2,...)
%
%   Input parameters:
%       x1,x2,...  - one or more vectors
%
%   Output parameters:
%       x1,x2,...  - input vectors as 3x1 vectors
%
%   POSITION_VECTOR(x1,x2,...) returns the given vectors x1,x2,... as
%   3x1 vectors. If the input is only two dimensional a zero is added as
%   the z dimension.
%
% see also: column_vector, row_vector
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ===================================
isargposition(varargin{:});


%% ===== Computation =====================================================
for ii = 1:nargin
    vec = row_vector(varargin{ii});
    if size(vec,2)<3
        varargout{ii} = [vec 0];
    else
        varargout{ii} = vec;
    end
end
