function isequallength(x1,varargin)
%ISEQUALSIZE tests if the given arrays have the same size
%   Usage: isequallength(x1,x2,...)
%
%   Input options:
%       x1,x2,...   - matrices to test
%
%   ISEQUALLENGTH(x1,x2,...) tests if all given arrays have the same size.
%
%   see also: isequalsize

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ===================================
if nargin<2
    error('%s: you need at least two input arguments.',upper(mfilename));
end


%% ===== Checking for equal size =========================================
for ii = 1:nargin-1
    if ~isequal(length(varargin{ii}),length(x1))
        error('%s and %s have not the same length',inputname(1),inputname(2));
    end
end
