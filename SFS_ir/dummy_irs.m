function irs = dummy_irs()
% DUMMY_IRS creates a dummy dirac pulse IR set
%
%   Usage: irs = dummy_irs()
%
%   Output parameters:
%       irs   - irs struct
%
%   DUMMY_IRS() creates a dummy IR data set (Dirac impulse) to check
%   processing without IRs. It has a resolution of 1 deg for phi and theta.
%
%   See also: new_irs, IR_format.txt

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


%% ===== Computation =====================================================
% create dirac pulse
nsamples = 1024;
ir = zeros(nsamples,1);
ir(300) = 1;
% angles of dummy irs
theta = rad(-90:89);
phi = rad(-180:179);
% replicate ir for all directions
ir = repmat(ir,1,length(phi)*length(theta));
% store data
irs = new_irs();
irs.left = ir;
irs.right = ir;
tmp = repmat(phi,length(theta),1);
irs.apparent_azimuth = tmp(:)';
irs.apparent_elevation = repmat(theta,1,length(phi));
irs.source_position = [0 2.333 0]';
irs.head_position = [0 0 0]';
irs.distance = 2.333;
irs.description = ['HRIR dummy set (Dirac pulse) for testing your',...
                   'frequency response, etc.'];
