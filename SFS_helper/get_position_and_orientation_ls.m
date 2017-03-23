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
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
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
