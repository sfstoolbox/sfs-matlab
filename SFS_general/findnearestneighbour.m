function [C,idx] = findnearestneighbour(A,b,number_of_neighbours)
%FINDNEARESTNEIGHBOUR finds the n nearest neighbours
%
%   Usage: [C,idx] = findnearestneighbour(A,b,[number_of_neighbours]);
%
%   Input parameters:
%       A                     - matrix
%       b                     - column to search for in A
%       number_of_neighbours  - number of nearest neighbours to find
%
%   output parameters:
%       C                     - found neighbour columns
%       idx                   - indices of found columns in matrix
%
%   FINDNEARESTNEIGHBOUR(A,b,number_of_neighbours) returns a number_of_neighbours
%   column vectors with the nearest neighbour points from the matrix A to the
%   point b. In addition to the values, the indices are also returned.
%
%   See also: find, findrows, get_ir

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Team                                   *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 3;
narginchk(nargmin,nargmax);
if nargin==nargmax-1
    number_of_neighbours = 1;
end
% Ensure column vector
if size(b,2)>1
    b=b';
end


%% ===== Computation =====================================================
% Calculate distance between points
distance = vector_norm(bsxfun(@minus,A,b),1);
% Sort the distances in order to find the n lowest once
[~,idx] = sort(distance);
idx = idx(1:min(number_of_neighbours,length(idx)));
C = A(:,idx);
