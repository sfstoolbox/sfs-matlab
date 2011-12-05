function [x,y,P] = wave_field_mono_line_source(X,Y,xs,f,conf)
%WAVE_FIELD_MONO_LINE_SOURCE simulates a wave field for a line source
%   Usage: [x,y,P] = wave_field_mono_line_source(X,Y,xs,f,conf)
%          [x,y,P] = wave_field_mono_line_source(X,Y,xs,f)
%
%   Input parameters:
%       X           - length of the X axis (m); single value or [xmin,xmax]
%       Y           - length of the Y axis (m); single value or [ymin,ymax]
%       xs          - position of line source (m)
%       f           - monochromatic frequency (Hz)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       x           - corresponding x axis
%       y           - corresponding y axis
%       P           - Simulated wave field
%
%   WAVE_FIELD_MONO_LINE_SOURCE(X,Y,xs,f,conf) simulates a wave 
%   field of a line source positioned at xs. 
%   To plot the result use plot_wavefield(x,y,P).
%
%   References:
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: plot_wavefield, wave_field_imp_point_source 

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 5;
error(nargchk(nargmin,nargmax,nargin));
isargvector(X,Y);
isargposition(xs);
xs = position_vector(xs);
isargpositivescalar(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
% Reference position for the amplitude (correct reproduction of amplitude
% at y = yref).
xref = conf.xref;
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
% Check if yref is in the given y space
% FIXME: write a function to check xref position
%if yref>max(y)
%    error('%s: yref has be smaller than max(y) = %.2f',...
%        upper(mfilename),max(y));
%end
% Create a x-y-grid to avoid a loop
[xx,yy] = meshgrid(x,y);
% Source model for a line source G_2D(x,omega)
P = line_source(xx,yy,xs,f);
% Scale signal (at yref)
P = norm_wave_field(P,x,y,conf);


% ===== Plotting =========================================================
if(useplot)
    plot_wavefield(x,y,P,conf);
end
