function varargout = axis_vector(varargin)
%AXIS_VECTOR returns 2x1 axis vectors [min max]
%
%   Usage: [x1,x2,...] = axis_vector(x1,x2,...)
%
%   Input parameters:
%       x1,x2,...  - one or more vectors
%
%   Output parameters:
%       x1,x2,...  - input vectors as 2x1 vectors
%
%   AXIS_VECTOR(x1,x2,...) returns the given vectors x1,x2,... as
%   2x1 vectors. If the input is only one dimensional, its value will be
%   replicated.
%
% see also: xyz_axes, xyz_grid

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


%% ===== Checking of input  parameters ===================================
isargvector(varargin{:});


%% ===== Computation =====================================================
for ii = 1:nargin
    vec = row_vector(varargin{ii});
    if size(vec,2)<2
        varargout{ii} = [vec vec];
    else
        varargout{ii} = vec;
    end
end
