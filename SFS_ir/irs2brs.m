function brs = irs2brs(irs)
%IRS2BRS converts a irs data set to an brs set suitable for the SSR
%
%   Usage: brs = irs2brs(irs)
%
%   Input parameters:
%       irs     - irs data set
%
%   Output parameters:
%       brs     - brs data set
%
%   IRS2BRS(irs) converts a irs data set into a brs set suitable for the
%   SoundScape Renderer. The brs data set is a matrix containing the
%   channels for all directions.
%
%   see also: set_wfs_25d

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
narginchk(nargmin,nargmax)
if nargin==nargmax-1
    conf = SFS_config;
end
check_irs(irs);
isargstruct(conf);


%% ===== Main ===========================================================

% Check if only one elevation angle is given
if length(unique(irs.apparent_elevation))~=1
    error(['%s: Your irs set has different elevation angles, which is',...
        ' not supported by the SoundScape Renderer.'],upper(mfilename));
end

% TODO: check the order of angles
%       I think the user have to check this by itself. Because the user
%       could also be interested in a particular angle order

for ii = 1:length(irs.apparent_azimuth)
    brs(:,(ii-1)*2+1:ii*2) = [irs.left(:,ii) irs.right(:,ii)];
end



