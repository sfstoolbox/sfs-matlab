function irs = dummy_irs()
% DUMMY_IRS creates a dummy dirac pulse IR set
%
%   Usage: irs = dummy_irs()
%
%   Output parameters:
%       irs   - irs struct
%
%   DUMMY_IRS() creates a dummy IR data set (Dirac impulse) to check
%   processing without IRs.
%
%   See also: new_irs, IR_format.txt

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


%% ===== Computation =====================================================
% Create dirac pulse
nsamples = 1024;
ir = zeros(nsamples,1);
ir(300) = 1;

irs = new_irs();
for ii=0:359
    irs.left(:,ii+1) = ir;
    irs.right(:,ii+1) = ir;
    irs.apparent_azimuth(ii+1) = correct_azimuth(ii/180*pi);
    irs.apparent_elevation(ii+1) = correct_elevation(0);
end
irs.description = ['HRIR dummy set (Dirac pulse) for testing your',...
                   'frequency response, etc.'];
% Reorder entries
irs = correct_irs_angle_order(irs);
