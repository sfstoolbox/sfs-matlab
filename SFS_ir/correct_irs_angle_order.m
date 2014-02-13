function irs = correct_irs_angle_order(irs,conf)
%CORRECT_IRS_ANGLE_ORDER reorders the angle entries of a irs to be increasing
%
%   Usage: irs = correct_irs_angle_order(irs,[conf])
%
%   Input options
%       irs     - irs struct
%       conf    - optional configuration struct (see SFS_config)
%
%   Output options
%       irs     - irs struct with corrected angle ordering
%
%   CORRECT_IRS_ANGLE_ORDER(irs,conf) corrects the order of the azimuth and
%   elevation entries to be increasing over the whole range. This is needed for
%   the interpolation functions.
%
%   See also: ir_intpol

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


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
if nargin==nargmax-1
    conf = SFS_config;
end
if conf.debug
    check_irs(irs);
end


%% ===== Computation =====================================================
% Sort azimuth angle and reorder the whole irs
irs.apparent_azimuth = correct_azimuth(irs.apparent_azimuth);
[~,idx] = sort(irs.apparent_azimuth);
irs.left = irs.left(:,idx);
irs.right = irs.right(:,idx);
irs.apparent_azimuth = irs.apparent_azimuth(idx);
irs.apparent_elevation = irs.apparent_elevation(idx);
if ~isequal(size(irs.head_azimuth),[1 1])
    irs.head_azimuth = irs.head_azimuth(idx);
end
if ~isequal(size(irs.head_elevation),[1 1])
    irs.head_elevation = irs.head_elevation(idx);
end
if ~isequal(size(irs.torso_azimuth),[1 1])
    irs.torso_azimuth = irs.torso_azimuth(idx);
end
if ~isequal(size(irs.torso_elevation),[1 1])
    irs.torso_elevation = irs.torso_elevation(idx);
end
if ~isequal(size(irs.distance),[1 1])
    irs.distance = irs.distance(idx);
end
if ~(isequal(size(irs.source_position),[3 1]) || ...
        isequal(size(irs.source_position),[1 3]))
    irs.source_position = irs.source_position(:,idx);
end
