function [x,y,P,ls_activity] = wave_field_mono_wfs_25d(X,Y,xs,ys,L,f,src,conf)
%WAVE_FIELD_MONO_WFS_25D simulates a wave field for 2.5D WFS
%   Usage: [x,y,P] = wave_field_mono_wfs_25d(X,Y,xs,ys,L,f,src,conf)
%          [x,y,P] = wave_field_mono_wfs_25d(X,Y,xs,ys,L,f,src)
%
%   Input parameters:
%       X           - length of the X axis (m); single value or [xmin,xmax]
%       Y           - length of the Y axis (m); single value or [ymin,ymax]
%       xs          - x position of point source (m)
%       ys          - y position of point source (m)
%       L           - array length (m)
%       f           - monochromatic frequency (Hz)
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs, ys are the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       x           - corresponding x axis
%       y           - corresponding y axis
%       P           - Simulated wave field
%
%   WAVE_FIELD_MONO_WFS_25D(X,Y,xs,ys,L,f,src,conf) simulates a wave 
%   field of the given source type (src) using a WFS 2.5 dimensional driving 
%   function in the temporal domain. This means by calculating the integral for 
%   P with a summation.
%   To plot the result use plot_wavefield(x,y,P).
%
%   References:
%       Spors2009 - Physical and Perceptual Properties of Focused Sources in
%           Wave Field Synthesis (AES127)
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: plot_wavefield, wf_SDM_25D

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 7;
nargmax = 8;
error(nargchk(nargmin,nargmax,nargin));
isargvector(X,Y);
isargscalar(xs,ys);
isargpositivescalar(L,f);
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


%% ===== Variables ======================================================

% Setting x- and y-axis
[X,Y] = setting_xy_ranges(X,Y,conf);

% Geometry
x = linspace(X(1),X(2),xysamples);
y = linspace(Y(1),Y(2),xysamples);


%% ===== Computation ====================================================

% Calculate the wave field in time-frequency domain
%
% Get the position of the loudspeakers and its activity
[x0,y0,phi] = secondary_source_positions(L,conf);
ls_activity = secondary_source_selection(x0,y0,phi,xs,ys,src);
% Generate tapering window
win = tapwin(L,ls_activity,conf);
ls_activity = ls_activity .* win;
% Create a x-y-grid to avoid a loop
[X,Y] = meshgrid(x,y);
% Initialize empty wave field
P = zeros(length(y),length(x));
% Integration over secondary source positions
for ii = 1:length(x0)

    % ====================================================================
    % Secondary source model G(x-x0,omega)
    % This is the model for the loudspeakers we apply. We use closed cabinet
    % loudspeakers and therefore point sources.
    G = point_source(X,Y,x0(ii),y0(ii),f);

    % ====================================================================
    % Driving function D(x0,omega)
    D = driving_function_mono_wfs_25d(x0(ii),y0(ii),phi(ii),xs,ys,f,src,conf);

    % ====================================================================
    % Integration
    %              /
    % P(x,omega) = | D(x0,omega) G(x-x0,omega) dx0
    %              /
    %
    % see: Spors2009, Williams1993 p. 36
    %
    % NOTE: win(ii) is the factor of the tapering window in order to have fewer
    % truncation artifacts. If you don't use a tapering window win(ii) will
    % always be one.
    P = P + win(ii)*D.*G;

end

% === Scale signal (at [xref yref]) ===
P = norm_wave_field(P,x,y,conf);


% ===== Plotting =========================================================
if(useplot)
    plot_wavefield(x,y,P,L,ls_activity,conf);
end
