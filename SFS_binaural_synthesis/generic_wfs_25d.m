function ir = generic_wfs_25d(xs,src,L,conf)
%GENRIC_WFS_25D Generate a IR for the generic renderer of the SSR
%
%   Usage: ir = generic_wfs_25d(xs,src,L,[conf])
%
%   Input parameters:
%       xs      - virtual source position [ys > Y0 => focused source] (m)
%       src     - source type: 'pw' -plane wave
%                              'ps' - point source
%                              'fs' - focused source
%       L       - Length of linear loudspeaker array (m)
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       ir      - Impulse response for the desired WFS array loudspeaker
%                 (nx1)
%
%   GENERIC_WFS_25D(xs,src,L,conf) calculates an impulse
%   response for a virtual source at xs for the loudspeakers of a WFS
%   array. every loudspeaker of the array is represented by one column in
%   the impulse response.
%
% see also: brs_wfs_25d, brs_point_source, auralize_ir

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

% FIXME: at the moment the first loudspeaker of your array has to be on the x-axis
% (which means phi=0). If you have another setup (like we have in Pinta) you
% have to manually edit the secondary_source_positions.m function in order to
% get the first loudspeaker at the desired location.


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
xs = position_vector(xs);
isargpositivescalar(L);
isargchar(src);

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
fs = conf.fs;                 % sampling frequency
c = conf.c;                   % speed of sound
N = conf.N;                   % target length of BRS impulse responses
useplot = conf.useplot;       % Plot results?


%% ===== Variables ======================================================
% Secondary sources
x0 = secondary_source_positions(L,conf);
x0 = secondary_source_selection(x0,xs,src);
nls = size(x0,1);
% generate tapering window
win = tapering_window(x0,conf);


%% ===== IR =============================================================

% === Dirac pulse ===
% Generate a dirac pulse to use in the delaying and add procedure
dirac = zeros(1024,1);
dirac(300) = 1;
lenir = length(dirac);
% Initial values
ir = zeros(N,nls);
dt = zeros(1,nls);
a = zeros(1,nls);

% Create a IR for every single loudspeaker
for n=1:nls

    % ====================================================================
    % Driving function to get weighting and delaying
    [a(n),delay] = ...
        driving_function_imp_wfs_25d(x0(n,:),xs,src,conf);
    % Time delay in samples for the given loudspeaker
    dt(n) = ceil( delay*fs ) + 500;

    % Check if the length of the IR (conf.N) is long enough for the
    % needed time delay dt(n)
    if N<(dt(n)+lenir)
        error(['%s: The length of the IR conf.N is not large enough ' ...
               'to handle the needed time delay dt(n).'],upper(mfilename));
    end

    % Generate impulse response for the desired loudspeaker
    ir(:,n) = [zeros(1,dt(n)) a(n)*win(n)*dirac' zeros(1,N-dt(n)-lenir)]';

end

%% ===== Pre-equalization ================================================
ir = wfs_preequalization(ir,conf);


%% ===== Plot WFS parameters =============================================
if(useplot)
    figure
    plot(x0(:,1),dt);
    title('delay (taps)');
    grid on;

    figure
    plot(x0(:,1),a);
    title('amplitude');
    grid on;
end

