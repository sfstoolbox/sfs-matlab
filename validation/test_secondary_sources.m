function boolean = test_secondary_sources(modus)
%TEST_SECONDARY_SOURCES tests the correctness of the secondary source
%functions
%
%   Usage: boolean = test_secondary_sources(modus)
%
%   Input parameters:
%       modus   - 0: numerical
%                 1: visual
%
%   Output parameters:
%       booelan - true or false
%
%   TEST_SECONDARY_SOURCES(modus) checks if the functions, that calculates
%   the secondary source positions and directions are working correctly.

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

% AUTHOR: Hagen Wierstorf
% $LastChangedDate: $
% $LastChangedRevision: $
% $LastChangedBy: $


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));


%% ===== Main ============================================================
conf = SFS_config;
L = 4;
conf.dx0 = 0.2;
boolean = true;
% reference values
x0_linear_ref = [
   -2.0000         0         0         0    1.0000         0
   -1.8000         0         0         0    1.0000         0
   -1.6000         0         0         0    1.0000         0
   -1.4000         0         0         0    1.0000         0
   -1.2000         0         0         0    1.0000         0
   -1.0000         0         0         0    1.0000         0
   -0.8000         0         0         0    1.0000         0
   -0.6000         0         0         0    1.0000         0
   -0.4000         0         0         0    1.0000         0
   -0.2000         0         0         0    1.0000         0
         0         0         0         0    1.0000         0
    0.2000         0         0         0    1.0000         0
    0.4000         0         0         0    1.0000         0
    0.6000         0         0         0    1.0000         0
    0.8000         0         0         0    1.0000         0
    1.0000         0         0         0    1.0000         0
    1.2000         0         0         0    1.0000         0
    1.4000         0         0         0    1.0000         0
    1.6000         0         0         0    1.0000         0
    1.8000         0         0         0    1.0000         0
    2.0000         0         0         0    1.0000         0
    ];
x0_circle_ref = [
    2.0054         0         0   -1.0000         0         0
    1.9954    0.1997         0   -0.9950   -0.0996         0
    1.9656    0.3974         0   -0.9802   -0.1981         0
    1.9163    0.5911         0   -0.9556   -0.2948         0
    1.8479    0.7789         0   -0.9215   -0.3884         0
    1.7611    0.9591         0   -0.8782   -0.4783         0
    1.6569    1.1297         0   -0.8262   -0.5633         0
    1.5362    1.2890         0   -0.7660   -0.6428         0
    1.4002    1.4356         0   -0.6982   -0.7159         0
    1.2503    1.5678         0   -0.6235   -0.7818         0
    1.0880    1.6845         0   -0.5425   -0.8400         0
    0.9149    1.7845         0   -0.4562   -0.8899         0
    0.7326    1.8667         0   -0.3653   -0.9309         0
    0.5431    1.9304         0   -0.2708   -0.9626         0
    0.3482    1.9749         0   -0.1736   -0.9848         0
    0.1499    1.9997         0   -0.0747   -0.9972         0
   -0.0500    2.0047         0    0.0249   -0.9997         0
   -0.2494    1.9898         0    0.1243   -0.9922         0
   -0.4462    1.9551         0    0.2225   -0.9749         0
   -0.6387    1.9009         0    0.3185   -0.9479         0
   -0.8248    1.8279         0    0.4113   -0.9115         0
   -1.0027    1.7367         0    0.5000   -0.8660         0
   -1.1706    1.6282         0    0.5837   -0.8119         0
   -1.3269    1.5036         0    0.6617   -0.7498         0
   -1.4700    1.3640         0    0.7331   -0.6802         0
   -1.5985    1.2108         0    0.7971   -0.6038         0
   -1.7111    1.0457         0    0.8533   -0.5214         0
   -1.8068    0.8701         0    0.9010   -0.4339         0
   -1.8844    0.6859         0    0.9397   -0.3420         0
   -1.9433    0.4948         0    0.9691   -0.2468         0
   -1.9830    0.2989         0    0.9888   -0.1490         0
   -2.0029    0.1000         0    0.9988   -0.0498         0
   -2.0029   -0.1000         0    0.9988    0.0498         0
   -1.9830   -0.2989         0    0.9888    0.1490         0
   -1.9433   -0.4948         0    0.9691    0.2468         0
   -1.8844   -0.6859         0    0.9397    0.3420         0
   -1.8068   -0.8701         0    0.9010    0.4339         0
   -1.7111   -1.0457         0    0.8533    0.5214         0
   -1.5985   -1.2108         0    0.7971    0.6038         0
   -1.4700   -1.3640         0    0.7331    0.6802         0
   -1.3269   -1.5036         0    0.6617    0.7498         0
   -1.1706   -1.6282         0    0.5837    0.8119         0
   -1.0027   -1.7367         0    0.5000    0.8660         0
   -0.8248   -1.8279         0    0.4113    0.9115         0
   -0.6387   -1.9009         0    0.3185    0.9479         0
   -0.4462   -1.9551         0    0.2225    0.9749         0
   -0.2494   -1.9898         0    0.1243    0.9922         0
   -0.0500   -2.0047         0    0.0249    0.9997         0
    0.1499   -1.9997         0   -0.0747    0.9972         0
    0.3482   -1.9749         0   -0.1736    0.9848         0
    0.5431   -1.9304         0   -0.2708    0.9626         0
    0.7326   -1.8667         0   -0.3653    0.9309         0
    0.9149   -1.7845         0   -0.4562    0.8899         0
    1.0880   -1.6845         0   -0.5425    0.8400         0
    1.2503   -1.5678         0   -0.6235    0.7818         0
    1.4002   -1.4356         0   -0.6982    0.7159         0
    1.5362   -1.2890         0   -0.7660    0.6428         0
    1.6569   -1.1297         0   -0.8262    0.5633         0
    1.7611   -0.9591         0   -0.8782    0.4783         0
    1.8479   -0.7789         0   -0.9215    0.3884         0
    1.9163   -0.5911         0   -0.9556    0.2948         0
    1.9656   -0.3974         0   -0.9802    0.1981         0
    1.9954   -0.1997         0   -0.9950    0.0996         0
    ];
x0_box_ref = [
   -2.0000    2.2000         0         0   -1.0000         0
   -1.8000    2.2000         0         0   -1.0000         0
   -1.6000    2.2000         0         0   -1.0000         0
   -1.4000    2.2000         0         0   -1.0000         0
   -1.2000    2.2000         0         0   -1.0000         0
   -1.0000    2.2000         0         0   -1.0000         0
   -0.8000    2.2000         0         0   -1.0000         0
   -0.6000    2.2000         0         0   -1.0000         0
   -0.4000    2.2000         0         0   -1.0000         0
   -0.2000    2.2000         0         0   -1.0000         0
         0    2.2000         0         0   -1.0000         0
    0.2000    2.2000         0         0   -1.0000         0
    0.4000    2.2000         0         0   -1.0000         0
    0.6000    2.2000         0         0   -1.0000         0
    0.8000    2.2000         0         0   -1.0000         0
    1.0000    2.2000         0         0   -1.0000         0
    1.2000    2.2000         0         0   -1.0000         0
    1.4000    2.2000         0         0   -1.0000         0
    1.6000    2.2000         0         0   -1.0000         0
    1.8000    2.2000         0         0   -1.0000         0
    2.0000    2.2000         0         0   -1.0000         0
    2.2000    2.0000         0   -1.0000         0         0
    2.2000    1.8000         0   -1.0000         0         0
    2.2000    1.6000         0   -1.0000         0         0
    2.2000    1.4000         0   -1.0000         0         0
    2.2000    1.2000         0   -1.0000         0         0
    2.2000    1.0000         0   -1.0000         0         0
    2.2000    0.8000         0   -1.0000         0         0
    2.2000    0.6000         0   -1.0000         0         0
    2.2000    0.4000         0   -1.0000         0         0
    2.2000    0.2000         0   -1.0000         0         0
    2.2000         0         0   -1.0000         0         0
    2.2000   -0.2000         0   -1.0000         0         0
    2.2000   -0.4000         0   -1.0000         0         0
    2.2000   -0.6000         0   -1.0000         0         0
    2.2000   -0.8000         0   -1.0000         0         0
    2.2000   -1.0000         0   -1.0000         0         0
    2.2000   -1.2000         0   -1.0000         0         0
    2.2000   -1.4000         0   -1.0000         0         0
    2.2000   -1.6000         0   -1.0000         0         0
    2.2000   -1.8000         0   -1.0000         0         0
    2.2000   -2.0000         0   -1.0000         0         0
    2.0000   -2.2000         0         0    1.0000         0
    1.8000   -2.2000         0         0    1.0000         0
    1.6000   -2.2000         0         0    1.0000         0
    1.4000   -2.2000         0         0    1.0000         0
    1.2000   -2.2000         0         0    1.0000         0
    1.0000   -2.2000         0         0    1.0000         0
    0.8000   -2.2000         0         0    1.0000         0
    0.6000   -2.2000         0         0    1.0000         0
    0.4000   -2.2000         0         0    1.0000         0
    0.2000   -2.2000         0         0    1.0000         0
         0   -2.2000         0         0    1.0000         0
   -0.2000   -2.2000         0         0    1.0000         0
   -0.4000   -2.2000         0         0    1.0000         0
   -0.6000   -2.2000         0         0    1.0000         0
   -0.8000   -2.2000         0         0    1.0000         0
   -1.0000   -2.2000         0         0    1.0000         0
   -1.2000   -2.2000         0         0    1.0000         0
   -1.4000   -2.2000         0         0    1.0000         0
   -1.6000   -2.2000         0         0    1.0000         0
   -1.8000   -2.2000         0         0    1.0000         0
   -2.0000   -2.2000         0         0    1.0000         0
   -2.2000   -2.0000         0    1.0000         0         0
   -2.2000   -1.8000         0    1.0000         0         0
   -2.2000   -1.6000         0    1.0000         0         0
   -2.2000   -1.4000         0    1.0000         0         0
   -2.2000   -1.2000         0    1.0000         0         0
   -2.2000   -1.0000         0    1.0000         0         0
   -2.2000   -0.8000         0    1.0000         0         0
   -2.2000   -0.6000         0    1.0000         0         0
   -2.2000   -0.4000         0    1.0000         0         0
   -2.2000   -0.2000         0    1.0000         0         0
   -2.2000         0         0    1.0000         0         0
   -2.2000    0.2000         0    1.0000         0         0
   -2.2000    0.4000         0    1.0000         0         0
   -2.2000    0.6000         0    1.0000         0         0
   -2.2000    0.8000         0    1.0000         0         0
   -2.2000    1.0000         0    1.0000         0         0
   -2.2000    1.2000         0    1.0000         0         0
   -2.2000    1.4000         0    1.0000         0         0
   -2.2000    1.6000         0    1.0000         0         0
   -2.2000    1.8000         0    1.0000         0         0
   -2.2000    2.0000         0    1.0000         0         0
   ];


% Calculate current values
% linear array
conf.array = 'linear';
x0_linear = secondary_source_positions(L,conf);
% circular array
conf.array = 'circle';
x0_circle = secondary_source_positions(L,conf);
% box form array
conf.array = 'box';
x0_box = secondary_source_positions(L,conf);

if modus==0
    % Numerical mode (quiet)
    if ~all(eq(size(x0_linear),size(x0_linear_ref))) || ...
            ~all(eq(size(x0_circle),size(x0_circle_ref))) || ...
            ~all(eq(size(x0_box),size(x0_box_ref))) || ...
            ~all(abs(x0_linear(:)-x0_linear_ref(:))<1e-4) || ...
            ~all(abs(x0_circle(:)-x0_circle_ref(:))<1e-4) || ...
            ~all(abs(x0_box(:)-x0_box_ref(:))<1e-4)
        boolean = false;
    end
elseif modus==1
    if ~all(eq(size(x0_linear),size(x0_linear_ref)))
        error('%s: wrong size of linear array.',upper(mfilename));
    elseif ~all(abs(x0_linear(:)-x0_linear_ref(:))<1e-4)
        error('%s: wrong value at linear array.',upper(mfilename));
    end
    if ~all(eq(size(x0_circle),size(x0_circle_ref)))
        error('%s: wrong size of circular array.',upper(mfilename));
    elseif ~all(abs(x0_circle(:)-x0_circle_ref(:))<1e-4)
        error('%s: wrong value at circular array.',upper(mfilename));
    end
    if ~all(eq(size(x0_box),size(x0_box_ref))) 
        error('%s: wrong size of box shaped array.',upper(mfilename));
    elseif ~all(abs(x0_box(:)-x0_box_ref(:))<1e-4)
        error('%s: wrong value at box shaped array.',upper(mfilename));
    end
elseif modus==2
    % Graphical mode
    close all;
    % ranges of the axis
    range = L/2+0.5;
    % draw results
    figure
    title('Linear loudspeaker array');
    axis([-1*range range -1*range range]);
    draw_loudspeakers(x0_linear,1,conf);
    figure
    title('Circular loudspeaker array');
    axis([-1*range range -1*range range]);
    draw_loudspeakers(x0_circle,1,conf);
    figure
    title('Box shape loudspeaker array');
    axis([-1*range range -1*range range]);
    draw_loudspeakers(x0_box,1,conf);
else
    error(['%s: modus has to be 0 (numerical quiet), 1 (numerical), ', ...
            'or 2 (graphical).'],upper(mfilename));
end
