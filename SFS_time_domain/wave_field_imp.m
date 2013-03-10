function [x,y,p] = wave_field_imp(X,Y,x0,d,t,conf)
%WAVE_FIELD_IMP returns the wave field in time domain of a loudspeaker array
%
%   Usage: [x,y,p] = wave_field_imp_nfchoa_25d(X,Y,x0,d,[conf])
%
%   Input options:
%       X           - length of the X axis (m); single value or [xmin,xmax]
%       Y           - length of the Y axis (m); single value or [ymin,ymax]
%       x0          - positions of secondary sources
%       d           - driving function of secondary sources
%       t           - time (samples)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output options:
%       x,y         - x- and y-axis of the wave field
%       p           - wave field (length(y) x length(x))
%
%   WAVE_FIELD_IMP(X,Y,x0,d,t,conf) computes the wave field synthesized by a 
%   loudspekaer array driven by individual driving functions to the time t.
%
%   To plot the result use:
%   conf.plot.usedb = 1;
%   plot_wavefield(x,y,p,conf);

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


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 6;
narginchk(nargmin,nargmax);
isargvector(X,Y);
isargsecondarysource(x0);
isargmatrix(d);
isargscalar(t);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end
if size(x0,1)~=size(d,2)
    error(['%s: The number of secondary sources (%i) and driving ', ...
        'signals (%i) does not correspond.'], ...
        upper(mfilename),size(x0,1),size(d,2));
end

%% ===== Configuration ==================================================
% Plotting result
useplot = conf.useplot;
% Speed of sound
c = conf.c;
% Sampling rate
fs = conf.fs;
% Debug mode
debug = conf.debug;


%% ===== Computation =====================================================

% Spatial grid
[xx,yy,x,y] = xy_grid(X,Y,conf);

% time reversal of driving function due to propagation of sound
% later parts of the driving function are emitted later by secondary
% sources
%d = d(end:-1:1,:);

% Apply bandbass filter
if(0)
    d=bandpass(d,conf);
end

% Initialize empty wave field
p = zeros(length(y),length(x));

% Integration over loudspeaker
for ii = 1:size(x0,1)

    % ================================================================
    % Secondary source model: Greens function g3D(x,t)
    % distance of secondary source to receiver position
    r = sqrt((xx-x0(ii,1)).^2 + (yy-x0(ii,2)).^2);
    % amplitude decay for a 3D monopole
    g = 1./(4*pi*r);

    % shift driving function
    %d(:,ii) = delayline(d(:,ii)',-size(d,1)+t,1,conf)';
    d(:,ii) = delayline(d(:,ii)',t,1,conf)';

    % Interpolate the driving function w.r.t. the propagation delay from
    % the secondary sources to a field point.
    % NOTE: the interpolation is required to account for the fractional
    % delay times from the loudspeakers to the field points
    t_vector = 1:length(d(:,ii));
    ds = interp1(t_vector,d(:,ii),r/c*fs,'spline');

    % ================================================================
    % Wave field p(x,t)
    p = p + ds .* g;

end

% === Checking of wave field ===
check_wave_field(p,t);


% === Plotting ===
if (useplot)
    conf.plot.usedb = 1;
    plot_wavefield(x,y,p,x0,conf);
end

% some debug stuff
if debug
    figure; imagesc(db(d)); title('driving functions'); caxis([-100 0]); colorbar;
    % figure; plot(win); title('tapering window');
    % figure; plot(delay*fs); title('delay (samples)');
    % figure; plot(weight); title('weight');
end
