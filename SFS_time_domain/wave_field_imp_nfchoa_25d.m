function [x,y,p] = wave_field_imp_nfchoa_25d(X,Y,xs,src,L,conf)
%WAVE_FIELD_IMP_NFCHOA_25D returns the wave field in time domain of an impulse
%
%   Usage: [x,y,p,ls_activity] = wave_field_imp_nfchoa_25d(X,Y,xs,src,L,[conf])
%
%   Input options:
%       X           - length of the X axis (m); single value or [xmin,xmax]
%       Y           - length of the Y axis (m); single value or [ymin,ymax]
%       xs          - position of point source (m)
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs, ys are the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%       L           - array length (m)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output options:
%       x,y         - x- and y-axis of the wave field
%       p           - wave field (length(y) x length(x))
%       ls_activity - activity of the secondary sources
%
%   WAVE_FIELD_IMP_NFCHOA_25D(X,Y,xs,src,L,conf) simulates a wave field of the
%   given source type (src) using a NFC-HOA 2.5 dimensional driving
%   function.
%   To plot the result use:
%   conf.plot.usedb = 1;
%   plot_wavefield(x,y,p,L,ls_activity,conf);

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


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 6;
narginchk(nargmin,nargmax);
isargvector(X,Y);
xs = position_vector(xs);
isargpositivescalar(L);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
% Plotting result
useplot = conf.useplot;
% Speed of sound
c = conf.c;
% Sampling rate
fs = conf.fs;
% Time frame to simulate
frame = conf.frame;
% Debug mode
debug = conf.debug;


%% ===== Computation =====================================================
% Get secondary sources
x0 = secondary_source_positions(L,conf);
nls = size(x0,1);

% Spatial grid
[xx,yy,x,y] = xy_grid(X,Y,conf);

% Calculate driving function
[d] = driving_function_imp_nfchoa_25d(x0,xs,src,L,conf);
% time reversal of driving function due to propagation of sound
% later parts of the driving function are emitted later by secondary
% sources
d = d(end:-1:1,:);

% shift driving function
for ii = 1:nls
    d(:,ii) = delayline(d(:,ii)',-size(d,1)+frame,1,conf)';
end

% Apply bandbass filter
if(0)
    d=bandpass(d,conf);
end

% Initialize empty wave field
p = zeros(length(y),length(x));

% Integration over loudspeaker
for ii = 1:nls

    % ================================================================
    % Secondary source model: Greens function g3D(x,t)
    % distance of secondary source to receiver position
    r = sqrt((xx-x0(ii,1)).^2 + (yy-x0(ii,2)).^2);
    % amplitude decay for a 3D monopole
    g = 1./(4*pi*r);

    % Interpolate the driving function w.r.t. the propagation delay from
    % the secondary sources to a field point.
    % NOTE: the interpolation is required to account for the fractional
    % delay times from the loudspeakers to the field points
    t = 1:length(d(:,ii));
    ds = interp1(t,d(:,ii),r/c*fs,'spline');

    % ================================================================
    % Wave field p(x,t)
    p = p + ds .* g;

end

% === Scale amplitude ===
% TODO: explain why
p = p ./ (pi*L);

% === Checking of wave field ===
check_wave_field(p,frame);


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
