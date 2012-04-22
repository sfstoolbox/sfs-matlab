function t = wave_front_time(X,xs,L,conf)
%WAVE_FRONT_TIME time of occurence of first echo for a linear WFS array
%
%   Usage: wave_front_time(X,xs,L,[conf])
%
%   Input parameters:
%       X       - listener position (m)
%       xs      - position of the virtual source (m)
%       L       - length of the linear loudspeaker array
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   WAVE_FRONT_TIME(X,xs,L) calculates the time of occuring of the first
%   echo (focused sources) resp. the last echo (virtual point source) at the
%   listener position X for a given point source location xs and a 
%   linear WFS loudspeaker array with a length of L.
%
%   see also: wave_front_direction, plot_wave_front_times

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


%% ===== Checking of input parameters ===================================
nargmin = 5;
nargmax = 6;
error(nargchk(nargmin,nargmax,nargin));
[X,xs] = position_vector(X,xs);
isargpositivescalar(L);

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
% Speed of sound
c = conf.c;


%% ===== Variables ======================================================
% Loudspeaker positions
x0 = secondary_source_positions(L,conf);
% Number of loudspeakers
nls = size(x0,1);


%% ===== Calculate a time axis ==========================================

% Geometry
%          x0(ii,:)                  X0
% x-axis <-^--^--^--^--^--^--^--^--^-|-^--^--^--^--^--^--^--^--^--
%             |                      |
%         R2 |  |                    |
%           |     | R                |      
%          x        |                |
%         xs          |              |
%                       O            |
%                       X            |
%                                    v
%                                  y-axis
%
% Distance between secondary source and listener: R
% Distance between secondary source and virtual source: R2

% Calculate arrival times at the virtual source position of the waves 
% emitted by the secondary sources
t2 = zeros(nls,1);
t = zeros(nls,1);
for i = 1:nls
    % Time between secondary sources and virtual source
    t2(ii) = norm(x0(ii,1:3)-xs)/c;
    % === Time, in which pre-echos occur ===
    t(ii) = norm(x0(ii,1:3)-X)/c;
end

% Adjust the time, so the virtual source arrives at time 0 at the listener
% position and use only the minimum time (focused sources).
t = min(t-1.*t2 - norm(xs-X)/c);
