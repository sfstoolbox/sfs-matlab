function [x,y,p,ls_activity] = wave_field_imp_wfs_25d(X,Y,xs,ys,L,src,conf)
%WAVE_FIELD_IMP_WFS_25D returns the wave field in time domain of an impulse
%
%   Usage: [x,y,p,ls_activity] = wave_field_imp_wfs_25d(X,Y,xs,ys,L,src,conf)
%          [x,y,p,ls_activity] = wave_field_imp_wfs_25d(X,Y,xs,ys,L,src)
%
%   Input options:
%       X           - length of the X axis (m); single value or [xmin,xmax]
%       Y           - length of the Y axis (m); single value or [ymin,ymax]
%       xs          - x position of point source (m)
%       ys          - y position of point source (m)
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
%   WAVE_FIELD_IMP_WFS_25D(X,Y,xs,ys,L,src,conf) simulates a wave 
%   field of the given source type (src) using a WFS 2.5 dimensional driving 
%   function with a delay line.
%   To plot the result use plot_wavefield(x,y,P).
%   NOTE: the pre-equalization filter is not integrated in this function at the
%   moment. But the main effect of this filter is to broaden the impulse a
%   little bit.

% AUTHOR: Hagen Wierstorf, Sascha Spors


%% ===== Checking of input  parameters ==================================
nargmin = 6;
nargmax = 7;
error(nargchk(nargmin,nargmax,nargin));
isargvector(X,Y);
isargscalar(xs,ys);
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

% Get loudspeaker positions, number of loudspeakers and their activity
[x0,y0,phi] = secondary_source_positions(L,conf);
nLS = length(x0);
ls_activity = secondary_source_selection(x0,y0,phi,xs,ys,src);
% Generate tapering window
win = tapwin(L,ls_activity,conf);
ls_activity = ls_activity .* win;

% Spatial grid
x = linspace(X(1),X(2),xysamples);
y = linspace(Y(1),Y(2),xysamples);
[X,Y] = meshgrid(x,y);

% Simulation for the given frame(s)
for sample = frame;
    % Initialize empty wave field
    p = zeros(length(y),length(x));
    % Integration over loudspeaker
    for l = 1:nLS

        % ================================================================
        % Driving function d2.5D(x0,t)
        [weight,delay] = driving_function_imp_wfs_25d(...
            x0(l),y0(l),phi(l),xs,ys,src,conf);

        % ================================================================
        % Secondary source model: Greens function g3D(x,t)
        % Distance of secondary source to receiver position
        r = sqrt((X-x0(l)).^2 + (Y-y0(l)).^2);
        % Greens function for a 3D monopole
        g = 1./(4*pi*r);

        % ================================================================
        % Driving function
        % NOTE: the interpolation is taking part because of the problem with the
        % sharp delta time points from the driving function and from the greens
        % function. Because we sample the time there is most often no overlap
        % between the two delta functions and the resulting wave field would be
        % zero.
        % FIXME: I should check, that I can combine the two delta functions in
        % the way I have done it here.
        % Calculate maximum time delay possible for the given axis size
        maxt = round(sqrt((X(end,end)-X(1,1))^2 + (Y(end,end)-Y(1,1))^2)/c*fs);
        % Add some additional offset
        maxt = maxt+500;
        % Create a time axis for the interpolation
        t = sample-maxt:sample+maxt;
        % create a driving function
        d = zeros(size(t));
        d(maxt+1) = weight * win(l);
        % Interpolate the driving function for the given delay time steps given
        % by the delta function from d, combined with the time steps given by
        % g.
        d = interp1(t,d,r/c*fs+delay*fs);

        % ================================================================
        % Wave field p(x,t)
        p = p + d .* g;

    end


    % === Checking of wave field ===
    check_wave_field(p);


    % === Plotting ===
    if (useplot)
        conf.plot.usedb = 1;
        plot_wavefield(x,y,p,L,ls_activity,conf);
    end

end
