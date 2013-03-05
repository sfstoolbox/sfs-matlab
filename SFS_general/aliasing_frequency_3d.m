function [fal] = aliasing_frequency_3d(x0,conf)
%ALIASING_FREQUENCY_3d returns the aliasing frequency for a 
% 3d spherical grid
%
%   Usage: fal = aliasing_frequency(x0,[conf])
%
%   Input options:
%       x0      - points of active secondary sources on a sphere
%       conf    - optional configuration struct (see SFS_config)
%
%   Output options:
%       fal     - aliasing frequency (Hz)
%
%   ALIASING_FREQUENCY(x0,conf) returns the aliasing frequency for the given
%   interspacing of secondary sources. The value is calculated.
%   
%   see also: wave_field_mono_wfs_3d

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

%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
if nargin==nargmax-1
    conf = SFS_config;
end
%[xs,X] = position_vector(xs,X);

%% ===== Configuration ==================================================
c = conf.c;

%% ===== Computation =====================================================
scalarproduct = (x0(1,:)*x0.')';
% calculate L2-norm
for ii = 1:length(x0)
    norm_vec(ii,1) = norm(x0(ii,:));
end
% calculate distance between points
cos_angle = scalarproduct./(norm_vec.*norm(x0(1,:)));
cos_angle = cos_angle(2:end);
dx0 = acos(max(cos_angle))*sqrt(x0(1,1)^2+x0(1,2)^2+x0(1,3)^2);
% calculate aliasing frequency
fal = c/(2*dx0);
