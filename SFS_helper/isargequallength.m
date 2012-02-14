function isargequallength(x1,varargin)
%ISARGEQUALSIZE tests if the given arrays have the same size
%   Usage: isargequallength(x1,x2,...)
%
%   Input options:
%       x1,x2,...   - matrices to test
%
%   ISARGEQUALLENGTH(x1,x2,...) tests if all given arigs have the same size.
%   Reports an error otherwise.
%
%   see also: isargequalsize

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


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
