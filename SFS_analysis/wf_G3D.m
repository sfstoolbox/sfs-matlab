function [x,y,G] = wf_G3D(X,Y,xs,ys,f,src,conf)
%WF_G3D simulates the wave field of a 3D Green's function 
%   Usage: [x,y,P] = wf_G3D(X,Y,xs,ys,f,src,conf)
%          [x,y,P] = wf_G3D(X,Y,xs,ys,f,src)
%
%   Input parameters:
%       X           - length of the X axis (m) [xaxis: -X/2:X/2]
%       Y           - length of the Y axis (m) [yaxis: -0.1:Y]
%       xs          - x position of point source (m)
%       ys          - y position of point source (m)
%       f           - monochromatic frequency (Hz)
%       src         - source type of the virtual source
%                         'pw' - plane wave
%                         'ps' - point source
%                         'fs' - point sink
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       x           - corresponding x axis
%       y           - corresponding y axis
%       G           - Simulated wave field
%
%   WF_G3D(X,Y,xs,ys,f,src,conf) simulates a wave field of the given 
%   source type (src) using the corresponding 3D Green's function in the 
%   temporal domain.  
%   To plot the result use plot_wavefield(x,y,G).
%
%   References:
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: plot_wavefield, wf_3DG_kx

% AUTHOR: Hagen Wierstorf


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
if ~isnumeric(f) || ~isscalar(f) || f<=0
    error('%s: f has to be a positive scalar!',upper(mfilename));
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
    conf = SFS_config;
end

% yref for amplitude scaling
yref = conf.yref;

% xy resolution
xysamples = conf.xysamples;

% phase of omega
phase = conf.phase;

% Plotting result
useplot = conf.useplot;

% Speed of sound
c = conf.c;


%% ===== Variables ======================================================

% General
omega = 2*pi*f;

% Geometry
y = linspace(-0.1,Y,xysamples);
x = linspace(-X/2,X/2,xysamples);


%% ===== Computation ====================================================

% Calculate the wave field in time-frequency domain 
%
% Create a x-y-grid to avoid a loop
[Y,X] = meshgrid(y,x);

if strcmp('ps',src)

    % === POINT SOURCE ===
    %
    % The 3D Green's function for a point source is given as
    % 
    %              1  e^(i w/c |x-xs|)
    % G(x-xs,w) = --- -----------------
    %             4pi      |x-xs|
    %
    % see: Williams1999, p. 198
    %
    % NOTE: the phase term e^(-i*phase) is only for the simulation of different 
    % time steps
    %
    G = 1/(4*pi) * exp(1i*omega/c.*sqrt((X-xs).^2+(Y-ys).^2)) ./ ...
            sqrt((X-xs).^2+(Y-ys).^2) .* exp(-1i*phase);


elseif strcmp('fs',src)

    % === POINT SINK ===
    %
    % The 3D Green's function for a point sink is given as
    % 
    %              1  e^(-i w/c |x-xs|)
    % G(x-xs,w) = --- ----------------
    %             4pi      |x-xs|
    %
    % see: Williams1999, p. 198
    %
    % NOTE: the phase term e^(-i*phase) is only for the simulation of different 
    % time steps
    %
    G = 1/(4*pi) * exp(-1i*omega/c.*sqrt((X-xs).^2+(Y-ys).^2)) ./ ...
            sqrt((X-xs).^2+(Y-ys).^2) .* exp(-1i*phase);


elseif strcmp('pw',src)

    % === PLANE WAVE ===
    %
    %FIXME: this works not correctly at the moment, why?
    %
    % The 3D green's function for a plane wave is given as
    %              
    % G = e^(i w/c nk x),
    %
    % where nk is a vector in the direction of k and x is the location vector.
    %
    % Use the position of the source as the direction vector for a plane
    % wave
    nxs = xs / sqrt(xs^2+ys^2);
    nys = ys / sqrt(xs^2+ys^2);
    %
    G = e^(1i*omega/c.*(nxs*X+nys*Y)) .* exp(-1i*phase);


else
    error('%s: the given sourcetype is not known!',upper(mfilename));
end

% === Scale signal (at xs,yref) ===
% Find index
[a,xidx] = find(x>xs,1);
[a,yidx] = find(y>yref,1);
% Scale signal to 1
G = 1*G/abs(G(xidx,yidx));


% ===== Plotting =========================================================
if(useplot)
    % Plot no loudspeakers
    conf.LSdist = 0.01;
    plot_wavefield(x,y,G,x(end),conf);
end
