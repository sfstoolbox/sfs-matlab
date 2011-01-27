function [x,y,P] = wf_SDM_25D(X,Y,xs,ys,L,f,src,conf)
%WF_SDM_25D simulates the wave field of a given source for 25D SDM
%   Usage: [x,y,P] = wf_SDM_25D(X,Y,xs,ys,L,f,src,conf)
%          [x,y,P] = wf_SDM_25D(X,Y,xs,ys,L,f,src)
%
%   Input parameters:
%       X           - length of the X axis (m) [xaxis: -X/2:X/2]
%       Y           - length of the Y axis (m) [yaxis: -0.1:Y]
%       xs          - x position of point source (m)
%       ys          - y position of point source (m)
%       L           - array length (m)
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
%       P           - Simulated wave field
%
%   WF_SDM_25D(X,Y,xs,ys,L,f,src,conf) simulates a wave field of the given
%   source type (src) using a SDM 2.5 dimensional driving function in the 
%   spectro-temporal domain. 
%   To plot the result use plot_wavefield(x,y,P).
%
%   NOTE: due to numerical problems with the fft and the bessel functions needed
%   in SDM (which resulted in an imaginary part which is hundreds of orders
%   greater/smaller than the real part) the FFT is done by hand in this
%   function. This results in a longer time to run this function. If you haven't
%   that time and you can life with the large argument approximation of the
%   bessel functions, which will result in a wrong evanescent part of the wave
%   field, you can use the wf_SDM_25D_approx function instead.
%
%   References:
%       Spors2010 - Reproduction of Focused Sources by the Spectral Division
%           Method
%       Spors2010 - Analysis and Improvement of Pre-equalization in
%       2.5-Dimensional Wave Field Synthesis
%
%   see also: plot_wavefield, wf_SDM_25D, fftx, ifftx

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 7;
nargmax = 8;
%error(nargchk(nargmin,nargmax,nargin));

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

% Array position (m)
X0 = conf.X0;
Y0 = conf.Y0;
% Loudspeaker distcane
% NOTE: if dLS <= dx, than we have a continous loudspeaker
dx = conf.LSdist;

% Method to calculate driving function (only for non-aliased part)
withev = conf.withev;  % with evanescent waves

% Reference position for the amplitude (correct reproduction of amplitude
% at y = yref).
yref = conf.yref;

% Use tapering window?
usetapwin = conf.usetapwin;

% xy resolution
xysamples = 300;

c = conf.c;
phase = conf.phase;

% Plotting
useplot = conf.useplot;


%% ===== Variables ======================================================

% General
omega = 2*pi*f;

% init variables
kxrep=2*pi/dx;
Nrep=6;             % number of spectral repetitions

% Aliasing condition
kxal = omega/c;
% Factor by which kx is extended of kx = omega/c criteria
Nkx=1.5;

%kx = linspace(-Nkx*kxal,Nkx*kxal,Nkx*2000);
kx = linspace(-Nkx*kxal,Nkx*kxal,Nkx*xysamples*10);
y  = linspace(-0.1,Y,xysamples);
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


%% ===== Spectrum of driving function =================================== 
Dkx = zeros(1,length(kx));

if strcmp('pw',src)

    % ===== PLANE WAVE ===================================================
    to_be_implemented(mfilename);

elseif strcmp('ps',src)

    % ===== POINT SOURCE =================================================
    to_be_implemented(mfilename);

elseif strcmp('fs',src)

    % ===== FOCUSED SOURCE ===============================================
    %
    % Exact driving function for 25D synthesis
    %
    % D_25D(kx,w) = e^(i kx xs) ...
    %                                   ____________
    %                         H0^(2)( \|(w/c)^2-kx^2 |yref-ys| )
    %                     / - --------------_-_-_-_-_-_---------, |kx|<|w/c|
    %                     |      H0^(2)( \|(w/c)^2-kx^2 yref ) 
    %                    <        ____________
    %                     | K0( \|kx^2-(w/c)^2 |yref-ys| )
    %                     \ ----------_-_-_-_-_-_---------,       |kx|>|w/c|
    %                          K0( \|kx^2-(w/c)^2 yref )
    %
    % NOTE: the phase term e^(i phase) is only there in order to be able
    %       to simulate different time steps
    %
    Dkx(idxpr) =  exp(1i*kx(idxpr)*xs) .* ...
        besselh(0,2,sqrt( (omega/c)^2 - kx(idxpr).^2 )*abs(yref-ys)) ./ ...
        besselh(0,2,sqrt( (omega/c)^2 - kx(idxpr).^2 )*abs(yref-Y0)) .* ...
        exp(1i*phase);
    if(withev)
        Dkx(idxev) =  exp(1i*kx(idxev)*xs) .* ...
            besselk(0,sqrt(kx(idxev).^2 - (omega/c).^2)*abs(yref-ys)) ./ ...
            besselk(0,sqrt(kx(idxev).^2 - (omega/c).^2)*abs(yref-Y0)) .* ...
            exp(1i*phase);
    end
    % Convolution with the truncation window
    % FIXME: this doesn't work with evanescent waves at the moment
    %w = L * sin(kx*L/(2))./(kx*L/(2));
    %Dkx = conv2(Dkx,w,'same');

else
    % No such source type for the driving function
    error('%s: src has to be one of "pw", "ps", "fs"!',upper(mfilename));
end


%% ===== Spectrum of secondary sources ================================== 
Gkx = zeros(length(kx),length(y));
% Green's function for a point source in the spectro-temporal domain (see
% Spors2010)
%                                  ____________
%                 / -i/4 H0^(2)( \|(w/c)^2-kx^2 y )
% G_3D(kx,y,w) = <                ____________
%                 \ 1/(2pi) K0( \|kx^2-(w/c)^2 y )
%
[K,Y] = meshgrid(kx(idxpr),abs(y-Y0));
Gkx(idxpr,:) = -1j/4 .* besselh(0,2,sqrt( (omega/c)^2 - K.^2 ).* Y)';
if(withev)
    [K,Y] = meshgrid(kx(idxev),abs(y-Y0));
    Gkx(idxev,:) = 1/(2*pi) .* besselk(0,sqrt( K.^2 - (omega/c)^2).* Y)';
end

%% ===== Reproduced wave field ========================================== 

% Pkx = Dkx * Gkx
Pkx = repmat(Dkx',1,length(y)) .* Gkx;


%% ===== Inverse spatial Fourier transformation =========================
% 
%            /
% P(x,y,w) = | Pkx(kx,y,w) * e^(-i kx x) dkx
%            /
%
P = zeros(length(x),length(y));
for n=1:length(x)
    % The following loop can be done faster by using the line below with repmat
    %for m=1:length(y)
    %    P(n,m) = sum ( Pkx(:,m) .* exp(-1j*kx*x(n))' );
    %end
    P(n,:) = sum ( Pkx .* repmat(exp(-1j*kx*x(n))',1,xysamples),1 );
end

% === Scale signal (at xs,yref) ===
% Find index
[a,xidx] = find(x>xs,1);
[a,yidx] = find(y>yref,1);
% Scale signal to 1
P = 1*P/abs(P(xidx,yidx));


%% ===== Plotting ========================================================
if(useplot)
    plot_wavefield(x,y,P,L,conf);
end
