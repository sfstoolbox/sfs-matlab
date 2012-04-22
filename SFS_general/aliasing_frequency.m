function fal = aliasing_frequency(dx0,conf)
%ALIASING_FREQUENCY returns the aliasing frequency
%
%   Usage: fal = aliasing_frequency(dx0,[conf])
%
%   Input options:
%       dx0     - distance between adjacent secondary sources (m)
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output options:
%       fal     - aliasing frequency (Hz)
%
%   ALIASING_FREQUENCY(dx0,conf) returns the aliasing frequency for the given
%   interspacing of secondary sources. The value is calculated after
%   spors2009.
%
%   S. Spors and J. Ahrens - Spatial sampling artifacts of wave field synthesis
%   for the reproduction of virtual point sources. 126th AES, May 2009.
%
%   see also: wave_field_mono_wfs_25d

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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox      sfs-toolbox@gmail.com *
%*****************************************************************************

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


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
% FIXME: better calculation possible?
fal = c/(2*dx0);
