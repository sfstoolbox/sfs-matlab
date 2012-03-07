function k = findrows(A,b,conf)
%FINDROWS finds indices of a given row within a matrix.
%
%   Usage: idx = findrows(A,b)
%
%   Input parameters:
%       A       - matrix
%       b       - vector to find in A
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameters:
%       idx     - indices of found columns in A
%
%   FINDROWS(A,b) returns a column vector with the indices of the rows
%   in the matrix A that are identical to the row vector b.  If no rows
%   in A are identical to b, an empty vector is returned.
%
%   The methods uses a for-loop, but it uses less memory and is in many
%   cases a lot faster than the vectorized methods
%
%      find( all( A == repmat(b, size(A, 1), 1), 2 ) )
%      find( all( A == b(ones(size(A, 1), 1),:), 2 ) )
%
%   See also: find, findcols
%

% AUTHOR: Peter John Acklam, Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));
if nargin==nargmax-1
    conf = SFS_config;
end
if conf.debug
    isargmatrix(A);
    isargvector(b);
end


%% ===== Computation =====================================================
k = find( A(:,1) == b(1) );
for j = 2:size(A, 2)
  k = k( A(k,j) == b(j) );
  if isempty(k)
     return
  end
end
