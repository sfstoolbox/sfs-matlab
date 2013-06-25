function [fal,dx0] = aliasing_frequency(x0,conf)
%ALIASING_FREQUENCY returns the aliasing frequency for the given secondary
%sources
%
%   Usage: [fal,dx0] = aliasing_frequency(x0,[conf])
%
%   Input options:
%       x0      - secondary sources / m
%       conf    - optional configuration struct (see SFS_config)
%
%   Output options:
%       fal     - aliasing frequency / Hz
%       dx0     - mean distance between secondary sources / m
%
%   ALIASING_FREQUENCY(x0,conf) returns the aliasing frequency for the given
%   secondary sources. First the mean distance dx0 between the secondary sources
%   is calculated, afterwards the aliasing frequency is calculated after Spors
%   (2009) as fal = c/(2*dx0).
%   For a calculation that includes the dependency on the listener position have
%   a look at Start (1997).
%   
%   References:
%       S. Spors and J. Ahrens - Spatial sampling artifacts of wave field
%       synthesis for the reproduction of virtual point sources. 126th AES,
%       May 2009.
%       E. Start - Direct Sound Enhancement by Wave Field Synthesis. TU Delft,
%       1997.
%
%   See also: wave_field_mono_wfs, secondary_source_positions

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
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
    % get distance between x01 and all secondary sources within x02
    dist = bsxfun(@minus,x02,x01);
    dist = vector_norm(dist,2);
    % get the smallest distance (which is the one to the next source)
    dx0(ii) = min(dist);
end
% get the mean distance between all secondary sources
dx0 = mean(dx0);
% calculate aliasing frequency
fal = c/(2*dx0);
