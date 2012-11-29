function t = wave_front_time(X,xs,src,L,conf)
%WAVE_FRONT_TIME time of occurence of wave fronts for a WFS array
%
%   Usage: wave_front_time(X,xs,src,L,[conf])
%
%   Input parameters:
%       X       - listener position (m)
%       xs      - position of the virtual source (m)
%       src     - source type of the virtual source
%                     'pw' - plane wave (xs, ys are the direction of the
%                            plane wave in this case)
%                     'ps' - point source
%                     'fs' - focused source
%       L       - length of the linear loudspeaker array
%       conf    - optional configuration struct (see SFS_config)
%
%   WAVE_FRONT_TIME(X,xs,src,L,conf) calculates the time of occurence of the
%   single wave fronts at te listener positions. The single wave fronts are 
%   omitted by the secondary sources.
%
%   see also: wave_front_direction, driving_function_imp_wfs_25d 

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


%% ===== Checking of input parameters ===================================
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
[X,xs] = position_vector(X,xs);
isargchar(src);
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

% Calculate arrival times at the listener position of the waves
% emitted by the secondary sources
t = zeros(nls,1);
delay = zeros(nls,1);
for ii = 1:nls
    % Time between secondary sources and listener
    t(ii) = norm(x0(ii,1:3)-X)/c;
    % Time delay of the secondary sources
    [weight,delay(ii)] = driving_function_imp_wfs_25d(x0(ii,:),xs,src,conf);
end

t = t+delay;
