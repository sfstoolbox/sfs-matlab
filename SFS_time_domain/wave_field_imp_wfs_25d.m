function [x,y,p,ls_activity] = wave_field_imp_wfs_25d(X,Y,xs,L,src,conf)
%WAVE_FIELD_IMP_WFS_25D returns the wave field in time domain of an impulse
%
%   Usage: [x,y,p,ls_activity] = wave_field_imp_wfs_25d(X,Y,xs,L,src,conf)
%          [x,y,p,ls_activity] = wave_field_imp_wfs_25d(X,Y,xs,L,src)
%
%   Input options:
%       X           - length of the X axis (m); single value or [xmin,xmax]
%       Y           - length of the Y axis (m); single value or [ymin,ymax]
%       xs          - position of point source (m)
%       L           - array length (m)
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs, ys are the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       conf        - optional configuration struct (see SFS_config)
%
%   Output options:
%       x,y         - x- and y-axis of the wave field
%       p           - wave field (length(y) x length(x))
%       ls_activity - activity of the secondary sources
%       
%   WAVE_FIELD_IMP_WFS_25D(X,Y,xs,L,src,conf) simulates a wave 
%   field of the given source type (src) using a WFS 2.5 dimensional driving 
%   function with a delay line.
%   To plot the result use plot_wavefield(x,y,P).

% AUTHOR: Hagen Wierstorf, Sascha Spors
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 6;
error(nargchk(nargmin,nargmax,nargin));
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
% xy resolution
xysamples = conf.xysamples;
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

% debug mode
debug=0;

%% ===== Computation =====================================================

% Check if virtual source is positioned at the right position.
%FIXME: this has to be done for all possible array forms!
%       Put it in an extra function
%check_geometry(xs,ys,L,X0,Y0,src);
%if strcmp('fs',src) && ys<=Y0
%    error('%s: ys has to be greater than Y0 for a focused source.', ...
%        upper(mfilename));
%elseif strcmp('ps',src) && ys>=Y0
%    error('%s: ys has to be smaller than Y0 for a point source.', ...
%        upper(mfilename));
%end

% Setting the x- and y-axis
[X,Y] = setting_xy_ranges(X,Y,conf);

% Get secondary sources
x0 = secondary_source_positions(L,conf);
ls_activity = secondary_source_selection(x0,xs,src);
% Generate tapering window
win = tapwin(L,ls_activity,conf);
ls_activity = ls_activity .* win;

% Spatial grid
x = linspace(X(1),X(2),xysamples);
y = linspace(Y(1),Y(2),xysamples);
[xx,yy] = meshgrid(x,y);

% Use only active secondary sources
x0 = x0(ls_activity>0,:);
nls = size(x0,1);
win = win(ls_activity>0);

% Calculate maximum time delay possible for the given axis size
maxt = round(sqrt((X(1)-X(2))^2+(Y(1)-Y(2))^2)/c*fs);
% Add some additional offset
aoffset=0;
maxt = maxt+aoffset;
% Create time axis for field interpolation
t = 0:maxt;

% Calculate pre-equalization filter if required
if(usehpre)
    hpre=wfs_prefilter(conf);
end

% Calculate driving function prototype
if(usehpre)
    d = [zeros(1,aoffset) hpre zeros(1,length(t)-length(hpre)-aoffset)];
else
    d = [zeros(1,aoffset) 1 zeros(1,length(t)-1-aoffset)];
end

% Apply bandbass filter to the prototype
%d=bandpass(d,conf);
%figure; freqz(d);
    
% In a first loop calculate the weight and delay values.
% This is done in an extra loop, because all delay values are needed to
% calculate the time frame to use for the wave field
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

if(debug)
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
    
    
    % save driving functions (debug)
    if(debug)
        dds(ii,:)=ds;
    end
    
    % Interpolate the driving function w.r.t. the propagation delay from
    % the secondary sources to a field point.
    % NOTE: the interpolation is taking part because of the problem with the
    % sharp delta time points from the driving function and from the greens
    % function. Because we sample in time there is most often no overlap
    % between the two delta functions and the resulting wave field would be
    % zero.
    ds = interp1(t,ds,r/c*fs,'cubic');
    %ds = interp1(t,ds,r/c*fs,'spline');
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
    plot_wavefield(x,y,p,L,ls_activity,conf);
end

% some debug stuff
if(debug)
    figure; imagesc(db(dds)); title('driving functions');
    figure; plot(win); title('tapering window');
    figure; plot(delay*fs); title('delay (samples)');
    figure; plot(weight); title('weight');
end