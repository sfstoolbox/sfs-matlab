function [x,y,p,ls_activity] = wave_field_time_domain(X,Y,xs,ys,L,src,conf)
%WAVE_FIELD_TIME_DOMAIN returns the wave field in time domain of an impulse
%
%   Usage: [x,y,p,ls_activity] = wave_field_time_domain(X,Y,xs,ys,L,src,conf)
%          [x,y,p,ls_activity] = wave_field_time_domain(X,Y,xs,ys,L,src)
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
%   WAVE_FIELD_TIME_DOMAIN(X,Y,xs,ys,L,src,conf) simulates a wave 
%   field of the given source type (src) using a WFS 2.5 dimensional driving 
%   function with a delay line.
%   To plot the result use plot_wavefield(x,y,P).

% AUTHOR: Hagen Wierstorf, Sascha Spors


%% ===== Checking of input  parameters ==================================
nargmin = 6;
nargmax = 7;
error(nargchk(nargmin,nargmax,nargin));

if ~isnumeric(X) || ~isvector(X)
    error('%s: X has to be a vector!',upper(mfilename));
end
if ~isnumeric(Y) || ~isvector(Y)
    error('%s: Y has to be a vector!',upper(mfilename));
end
if ~isnumeric(xs) || ~isscalar(xs)
    error('%s: xs has to be a scalar!',upper(mfilename));
end
if ~isnumeric(ys) || ~isscalar(ys)
    error('%s: ys has to be a scalar!',upper(mfilename));
end
if ~isnumeric(L) || ~isscalar(L) || L<=0
    error('%s: L has to be a positive scalar!',upper(mfilename));
end
if ~ischar(src)
    error('%s: src has to be a string!',upper(mfilename));
end
if nargin<nargmax
    useconfig = true;
elseif ~isstruct(conf)
    error('%s: conf has to be a struct.',upper(mfilename));
else
    useconfig = false;
end


%% ===== Configuration ==================================================

% Load default configuration values
if(useconfig)
    disp('yes');
    conf = SFS_config;
end

% Array position (m)
X0 = conf.X0;
Y0 = conf.Y0;
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

% Loudspeaker distance
% NOTE: if dLS <= dx, than we have a continous loudspeaker
dLS = conf.LSdist;
% Reference position for the amplitude (correct reproduction of amplitude
% at y = yref).
yref = conf.yref;
% Use tapering window?
usetapwin = conf.usetapwin;
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
% Create movie
%usemovie = conf.usemovie;


%% ===== Computation =====================================================

% Setting the x- and y-axis
[X,Y] = setting_xy_ranges(X,Y,conf);

% Get loudspeaker positions for the given array
[x0,y0,phi] = secondary_source_positions(L,conf);
nLS = length(x0);

% Driving function and secondary source activity
[d,ls_activity] = driving_function_wfs(X,Y,xs,ys,L,src,conf);

% Spatial grid
x = linspace(X(1),X(2),xysamples);
y = linspace(Y(1),Y(2),xysamples);
[X,Y] = meshgrid(x,y);


for t0 = frame;

    p = zeros(length(y),length(x));
    d_reshaped = zeros(length(y),length(x));

    for l = 1 : nLS

        % distance of secondary source to receiver position
        r = sqrt( (X-x0(l)).^2 + (Y-y0(l)).^2 );
        % propagation duration from secondary source to receiver position
        t = (t0/fs - r./c) .* fs + 1; % in samples

        % shift d appropriately and interpolate between samples to obtain
        % the amplitude at the receiver position
        d_reshaped = reshape( interp1( (1:size(d,1)), d(:,l), t, 'linear'), ...
                                           xysamples, xysamples );

        % add to reproduced sound field (and assure causality)
        idx = ((r>0));
        p(idx) = p(idx) + d_reshaped(idx) ./ r(idx);

    end


    % === Plotting ===
    if (useplot)

        plot_wavefield(x,y,p,L,ls_activity,conf);

    end

    %if ( movie )
    %    GraphDefaults('two_talk');
    %    file_name = sprintf('plane_wave/%0.3d.png', frame);
    %    print('-dpng', file_name);
    %    frame = frame + 1;
    %else
    %    GraphDefaults('paper');
    %end

end


