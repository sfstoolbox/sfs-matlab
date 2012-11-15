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
