function [fal,dx0] = aliasing_frequency_3d(x0,conf)
%ALIASING_FREQUENCY_3d returns the aliasing frequency for a 
% 3d spherical grid
%
%   Usage: fal = aliasing_frequency(x0,[conf])
%
%   Input options:
%       x0      - secondary sources / m
%       conf    - optional configuration struct (see SFS_config)
%
%   Output options:
%       fal     - aliasing frequency / Hz
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

%% ===== Configuration ==================================================
c = conf.c;

%% ===== Computation =====================================================
% calculate the distance to the nearest secondary source for all secondary
% sources
for ii=1:size(x0,1)
    % first secondary source position
    x01 = x0(ii,1:3);
    % all other positions
    x02 = [x0(1:ii-1,1:3); x0(ii+1:end,1:3)];
    % get distance between points
    dist = bsxfun(@minus,x02,x01);
    dist = vector_norm(dist,2);
    % get the smallest distance
    dx0(ii) = min(dist);
end
% get the mean distance between all speakers
dx0 = mean(dx0);
% calculate aliasing frequency
fal = c/(2*dx0);
