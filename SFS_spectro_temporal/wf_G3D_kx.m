function [x,y,G] = wf_G3D_kx(X,Y,xs,ys,f,src,conf)
%WF_G3D_KX simulates the wave field of a 3D green's function using the kx space
%   Usage: [x,y,G] = wf_G3D_kx(X,Y,xs,ys,f,src,conf)
%          [x,y,G] = wf_G3D_kx(X,Y,xs,ys,f,src)
%
%   Input parameters:
%       X           - length of the X axis (m) [xaxis: -X/2:X/2]
%       Y           - length of the Y axis (m) [yaxis: -0.1:Y]
%       xs          - x position of point source (m)
%       ys          - y position of point source (m)
%       f           - Monochromatic frequency (Hz)
%       src         - source type of the virtual source
%                         'pw' - plane wave
%                         'ps' - point source
%                         'fs' - focused source
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       x           - corresponding x axis
%       y           - corresponding y axis
%       G           - Simulated wave field
%
%   WF_G3D_KX(X,Y,xs,ys,f,src,conf) simulates a wave field of the given
%   source type (src) using ithe corresponding Green's function in the 
%   spectro-temporal domain. 
%   To plot the result use plot_wavefield(x,y,G).
%
%   References:
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: plot_wavefield, wf_G3D

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
    error('%s: sourcetype has to be a string!',upper(mfilename));
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

% Method to calculate driving function (only for non-aliased part)
withev = conf.withev;  % with evanescent waves

% Reference position for the amplitude (correct reproduction of amplitude
% at y = yref).
yref = conf.yref;

% xy resolution
xysamples = conf.xysamples;

c = conf.c;
phase = conf.phase;

% Plotting
useplot = conf.useplot;


%% ===== Variables ======================================================

% General
omega = 2*pi*f;

% Aliasing condition
kxal = omega/c;
% Factor by which kx is extended of kx = omega/c criteria
Nkx=1.5;

%kx = linspace(-Nkx*kxal,Nkx*kxal,Nkx*10000);
kx = linspace(-Nkx*kxal,Nkx*kxal,Nkx*xysamples);
y  = linspace(0,Y,xysamples);
x  = linspace(-X/2,X/2,xysamples);

% Indexes for evanescent contributions and propagating part of the wave field
idxpr = (( abs(kx) <= (omega/c) ));
idxev = (( abs(kx) > (omega/c) ));


%% ===== Computation =====================================================

% Check if yref is in the given y space
if yref>max(y)
    error('%s: yref has be smaller than max(y) = %.2f',...
        upper(mfilename),max(y));
end

Gkx = zeros(length(kx),length(y));

if strcmp('ps',src)

    % === POINT SOURCE ===
    % Green's function for a point source in the spectro-temporal domain (see
    % Spors2010)
    %                                  ____________
    %                 / -i/4 H0^(1)( \|(w/c)^2-kx^2 y )
    % G_3D(kx,y,w) = <                ____________
    %                 \ 1/(2pi) K0( \|kx^2-(w/c)^2 y )
    %
    [K,Y] = meshgrid(kx(idxpr),abs(y-ys));
    Gkx(idxpr,:) = -1j/4 .* besselh(0,1,sqrt( (omega/c)^2 - K.^2 ).* Y)';
    if(withev)
        [K,Y] = meshgrid(kx(idxev),abs(y-ys));
        Gkx(idxev,:) = 1/(2*pi) .* besselk(0,sqrt( K.^2 - (omega/c)^2).* Y)';
    end


elseif strcmp('fs',src)

    % === POINT SINK ===
    % Green's function for a point sink in the spectro-temporal domain (see
    % Spors2010)
    %                                  ____________
    %                 / -i/4 H0^(2)( \|(w/c)^2-kx^2 y )
    % G_3D(kx,y,w) = <                ____________
    %                 \ 1/(2pi) K0( \|kx^2-(w/c)^2 y )
    %
    [K,Y] = meshgrid(kx(idxpr),abs(y-ys));
    Gkx(idxpr,:) = -1j/4 .* besselh(0,2,sqrt( (omega/c)^2 - K.^2 ).* Y)';
    if(withev)
        [K,Y] = meshgrid(kx(idxev),abs(y-ys));
        Gkx(idxev,:) = 1/(2*pi) .* besselk(0,sqrt( K.^2 - (omega/c)^2).* Y)';
    end


elseif strcmp('pw',src)

    % === PLANE WAVE ===
    % Green's function for a plane wave in the spectro-temporal domain (see
    % Ahrens2010)
    %
    % G_3D(kx,y,w) = 4 pi^2 delta(kx - w/c nxs) e^(i w/c nys y)
    %
    % Use the position of the source as the direction vector for a plane wave
    nxs = xs / sqrt(xs^2+ys^2);
    nys = ys / sqrt(xs^2+ys^2);
    % Simulate the delta function by only calculating that index, that is
    % different from zero
    idx = find(kx>omega/c*nxs,1);
    Gkx(idx,:) = 4*pi^2 .* exp(1i*omega/c*nys*y) .* exp(-1i*phase);

else
    error('%s: the given sourcetype is not known!',upper(mfilename));
end


%% ===== Inverse spatial Fourier transformation =========================
% 
%            /
% G(x,y,w) = | Gkx(kx,y,w) * e^(-i kx x) dkx
%            /
%
G = zeros(length(x),length(y));
for n=1:length(x)
    for m=1:length(y)
        G(n,m) = sum ( Gkx(:,m) .* exp(-1j*kx*x(n))' );
    end
end

% === Scale signal (at xs,yref) ===
% Find index
[a,xidx] = find(x>xs,1);
[a,yidx] = find(y>yref,1);
% Scale signal to 1
G = 1*G./abs(G(xidx,yidx));


%% ===== Plotting ========================================================
if(useplot)
    plot_wavefield(x,y,G,conf);
end
