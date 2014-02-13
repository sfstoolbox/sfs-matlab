function irspart = slice_irs(irs,idx)
%SLICE_IRS returns a part of an IRs set given by idx
%
%   Usage: irspart = slice_irs(irs,idx)
%
%   Input parameters:
%       irs     - IR data set
%       idx     - idx to slice out of the IR set
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       irspart - HRIR/BRIR containing only the part of the original IR set
%                 given by idx
%
%   SLICE_IRS(irs,idx) returns a part of the IR set irs given by idx. The new
%   part is a full IR set containing all the necessary struct elements.
%
%   see also: get_ir, read_irs

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
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);


%% ===== Slicing the IR set ==============================================
irspart = irs;
irspart.left = irs.left(:,idx);
irspart.right = irs.right(:,idx);
irspart.apparent_azimuth = irs.apparent_azimuth(idx);
irspart.apparent_elevation = irs.apparent_elevation(idx);
if ~all(size(irs.head_azimuth)==[1 1])
    irspart.head_azimuth = irs.head_azimuth(idx);
end
if ~all(size(irs.head_elevation)==[1 1])
    irspart.head_elevation = irs.head_elevation(idx);
end
if ~all(size(irs.torso_azimuth)==[1 1])
    irspart.torso_azimuth = irs.torso_azimuth(idx);
end
if ~all(size(irs.torso_elevation)==[1 1])
    irspart.torso_elevation = irs.torso_elevation(idx);
end
if ~all(size(irs.distance)==[1 1])
    irspart.distance = irs.distance(idx);
end
if ~all(size(irs.source_position)==[3 1])
    irspart.source_position = irs.source_position(idx);
end
if ~all(size(irs.source_reference)==[3 1])
    irspart.source_reference = irs.source_reference(idx);
end
