function [x,y,p,win,dds] = wave_field_imp_wfs_25d(X,Y,xs,src,L,conf)
%WAVE_FIELD_IMP_WFS_25D returns the wave field in time domain of an impulse
%
%   Usage: [x,y,p,ls_activity] = wave_field_imp_wfs_25d(X,Y,xs,src,L,[conf])
%
%   Input options:
%       X           - length of the X axis (m); single value or [xmin,xmax]
%       Y           - length of the Y axis (m); single value or [ymin,ymax]
%       xs          - position of point source (m)
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs, ys are the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       L           - array length (m)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output options:
%       x,y         - x- and y-axis of the wave field
%       p           - wave field (length(y) x length(x))
%       ls_activity - activity of the secondary sources
%
%   WAVE_FIELD_IMP_WFS_25D(X,Y,xs,src,L,conf) simulates a wave field of the
%   given source type (src) using a WFS 2.5 dimensional driving function with
%   a delay line.
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
% Use pre-equalization filter
usehpre = conf.usehpre;
% Debug mode
debug = conf.debug;
xref = position_vector(conf.xref);


%% ===== Computation =====================================================
% Get secondary sources
x0 = secondary_source_positions(L,conf);
x0 = secondary_source_selection(x0,xs,src,xref);
nls = size(x0,1);
% Generate tapering window
win = tapering_window(x0,conf);

% Spatial grid
[xx,yy,x,y] = xy_grid(X,Y,conf);

% Calculate maximum time delay possible for the given axis size
maxt = round(sqrt((X(1)-X(2))^2+(Y(1)-Y(2))^2)/c*fs);
% Add some additional pre-offset
aoffset=128;
maxt = maxt+aoffset;
% Create time axis for field interpolation
t = 0:maxt;

% Calculate pre-equalization filter if required
if usehpre
    hpre = wfs_prefilter(conf);
else
    hpre = 1;
end
% Calculate driving function prototype
d = [zeros(1,aoffset) hpre zeros(1,length(t)-length(hpre)-aoffset)];

% Apply bandbass filter to the prototype
if(0)
    d=bandpass(d,conf);
end

% In a first loop calculate the weight and delay values.
% This is done in an extra loop, because all delay values are needed to
% calculate the time frame to use for the wave field
% FIXME: the calculation of the right frame is not working correctly at the
% moment.
delay = zeros(nls,1);
weight = zeros(nls,1);
for ii = 1:nls
    % ================================================================
    % Driving function d2.5D(x0,t)
    [weight(ii),delay(ii)] = driving_function_imp_wfs_25d(x0(ii,:),xs,src,conf);
end
dmin=min(delay);

% If no explizit time frame is given calculate one
if isempty(frame)
    % Use only those delays for the calculation, that correspond to secondary
    % sources within the shown listening area
    idx = abs(x0(:,1))<max(abs(X(:))) & abs(x0(:,2))<max(abs(Y(:)));
    % If we haven't found any idx, use all entries
    if isempty(idx)
        idx = ones(nls,1);
    end
    % Get frame
    frame = max(round(delay(idx)*fs)) + 100;
end


% In a second loop simulate the wave field
% Initialize empty wave field
p = zeros(length(y),length(x));

if debug
    dds = zeros(nls,length(d));
end

% Integration over loudspeaker
for ii = 1:nls

    % ================================================================
    % Secondary source model: Greens function g3D(x,t)
    % distance of secondary source to receiver position
    r = sqrt((xx-x0(ii,1)).^2 + (yy-x0(ii,2)).^2);
    % amplitude decay for a 3D monopole
    g = 1./(4*pi*r);

    % ================================================================
    % Shift and weight prototype driving function
    % - less delay in driving function is more propagation time in sound
    %   field, hence the sign of the delay has to be reversed in the
    %   argument of the delayline function
    % - the proagation time from the source to the nearest secondary source
    %   is removed
    % - the main pulse in the driving function is shifted by aoffset and
    %   frame
    ds = delayline(d,frame-(delay(ii)-dmin)*fs,weight(ii)*win(ii),conf);

    % remember driving functions (debug)
    if debug
        dds(ii,1:length(ds)) = ds;
    end

    % Interpolate the driving function w.r.t. the propagation delay from
    % the secondary sources to a field point.
    % NOTE: the interpolation is required to account for the frcational
    % delay times from the loudspeakers to the field points
    t = 1:length(ds);
    ds = interp1(t,ds,r/c*fs,'spline');
    %ds = interp1(t,ds,r/c*fs,'cubic');
    %ds = interp1(t,ds,r/c*fs,'linear');
    %ds = interp1(t,ds,r/c*fs,'nearest');

    % ================================================================
    % Wave field p(x,t)
    p = p + ds .* g;

end

% === Checking of wave field ===
check_wave_field(p,frame);

% === Plotting ===
if (useplot)
    conf.plot.usedb = 1;
    plot_wavefield(x,y,p,x0,win,conf);
end

% some debug stuff
if debug
    figure; imagesc(db(dds)); title('driving functions'); caxis([-100 0]); colorbar;
    % figure; plot(win); title('tapering window');
    % figure; plot(delay*fs); title('delay (samples)');
    % figure; plot(weight); title('weight');
end
