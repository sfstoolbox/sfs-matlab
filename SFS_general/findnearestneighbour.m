function [C,idx] = findnearestneighbour(A,b,number_of_neighbours)
%FINDNEARESTNEIGHBOUR finds the n nearest neighbours
%
%   Usage: [C,idx] = findnearestneighbour(A,b,[number_of_neighbours]);
%
%   Input parameters:
%       A                     - matrix
%       b                     - colum to search for in A
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
%   See also: find, findrows

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
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
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 3;
narginchk(nargmin,nargmax);
if nargin==nargmax-1
    number_of_neighbours = 1;
end
% ensure column vector
if size(b,2)>1
    b=b';
end
if number_of_neighbours>size(A,2)
    error(['%s: your number of neighbours is larger than the available ', ...
        'points.'],upper(mfilename));
end


%% ===== Computation =====================================================
% calculate distance between points
distance = vector_norm(bsxfun(@minus,A,b),1);
% sort the distances in order to find the n lowest once
[~,idx] = sort(distance);
idx = idx(1:number_of_neighbours);
C = A(:,idx);
