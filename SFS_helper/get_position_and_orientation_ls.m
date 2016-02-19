function [xs,nxs] = get_position_and_orientation_ls(xs,conf);
%GET_POSITION_AND_ORIENTATION_LS returns [nx3] position and [nx3] orientation
%for a line source
%
%   Usage: [xs,nxs] = get_position_and_orientation_ls(xs,conf);
%
%   Input options:
%       xs          - combined position and orientation / m [nx3] or [nx6]
%       args        - list of args
%
%   Output options:
%       xs          - position of line source / m [nx3]
%       nxs         - orientation of line source / m [nx3]
%
%   GET_POSITION_AND_ORIENTATION_LS(xs,conf) returns the position and
%   orientation of a point source. The orientation is a vector perpendicular to
%   the traveling direction of the line source. For 2D or 2.5D the orientation
%   will always be returned as [0 0 1]. For 3D [0 0 1] will be returned if no
%   explicit orientation is given.
%
%   See also: secondary_source_selection, driving_function_mono_wfs

%*****************************************************************************
% Copyright (c) 2010-2016 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2016 Institut fuer Nachrichtentechnik                   *
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


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargnumeric(xs);
isargstruct(conf);


%% ===== Configuration ===================================================
dimension = conf.dimension;


%% ===== Checking for vector =============================================
% Handling of line source orientation
if (strcmp('2D',dimension) ||  strcmp('2.5D',dimension))
    % Ignore orientation for 2D and 2.5D
    if size(xs,2)>3
        warning('%s: %s-WFS ignores virtual line source orientation.', ...
            upper(mfilename),dimension);
    end
    xs(:,3) = 0;
    nxs = repmat([0 0 1],[size(xs,1) 1]);
else
    if size(xs,2)~=6
        warning('%s: set ls orientation to [0 0 1]',upper(mfilename));
        nxs = repmat([0 0 1],[size(xs,1) 1]);
    else
        nxs = xs(:,4:6) / norm(xs(1,4:6),2);
    end
end
xs = xs(:,1:3);
