function [C,idx] = findnearestneighbour(A,b,number_of_neighbours)
%FINDCOLS finds indices of a given column within a matrix.
%
%   Usage: idx = findcols(A,b)
%
%   Input parameters:
%       A       - matrix
%       b       - colum to search for in A
%
%   output parameters:
%       idx     - indices of found columns in matrix
%
%   FINDCOLS(A,b) returns a row vector with the indices of the columns
%   in the matrix A that are identical to the column vector b.  If no
%   columns in A are identical to b, an empty vector is returned.
%
%   See also: find, findrows

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************

% AUTHOR: Peter John Acklam


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 3;
narginchk(nargmin,nargmax);
if nargin==nargmax-1
    number_of_neighbours = 1;
end
b = column_vector(b);
if number_of_neighbours>size(A,2)
    error(['%s: your number of neighbours is larger than the available ', ...
        'points.'],upper(mfilename));
end


%% ===== Computation =====================================================
% loop through the points
distance = zeros(1,size(A,2));
for ii=1:size(A,2)
    distance(ii) = norm(A(:,ii)-b);
end
% sort the distances in order to find the n lowest once
[~,idx] = sort(distance);
idx = idx(1:number_of_neighbours);
C = A(:,idx);
