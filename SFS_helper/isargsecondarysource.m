function isargsecondarysource(varargin)
%ISARGSECONDARYSOURCE tests if the given arg is a matrix containing secondary
%   source positions and directions and return an error otherwise
%
%   Usage: isargsecondarysource(args)
%
%   Input options:
%       args        - list of args
%
%   ISARGSECONDARYSOURCE(args) tests if all given args are a matrix
%   containing secondary source positions and directions and returns
%   an error otherwise.
%
%   see also: isargmatrix

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


%% ===== Checking for vector =============================================
for ii = 1:nargin
    x0 = varargin{ii};
    if ~isnumeric(x0) || ndims(x0)~=2 || size(x0,2)~=6
        error(['%s need to be a nx6 matrix containing the secondary ', ...
            'sources positions and directions.'],inputname(ii));
    end
    for jj=1:size(x0,1)
        if abs(norm(x0(jj,4:6))-1)>1e-10
            error(['The norm of the direction of %s is not 1. ', ...
                'You can use the direction_vector function to get ', ...
                'correct direction values.'],inputname(ii));
        end
    end
end
