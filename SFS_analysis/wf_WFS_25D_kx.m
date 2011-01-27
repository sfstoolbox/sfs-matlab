function [x,y,P] = wf_WFS_25D_kx(X,Y,xs,ys,L,f,src,conf)
%WF_WFS_25D_KX simulates the wavefield for the 2.5D case using WFS
%   Usage: [x,y,P] = wf_WFS_25D_kx(X,Y,xs,ys,L,f,src,conf)
%          [x,y,P] = wf_WFS_25D_kx(X,Y,xs,ys,L,f,src)
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
%   wf_WFS_25D_kx(X,Y,xs,ys,L,f,src,conf) simulates a wave field of a virtual
%   source with the given source type (src) using a WFS 2.5 dimensional
%   driving function in the spectro-temporal domain.
%   To plot the result use plot_wavefield(x,y,P).
%
%   References:
%       Spors2009 - Physical and Perceptual Properties of Focused Sources in
%           Wave Field Synthesis
%       Spors2010 - Analysis and Improvement of Pre-equalization in
%       2.5-Dimensional Wave Field Synthesis
%
%   see also: plot_wavefield, wf_SDM_25D_kx, fftx, ifftx

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 7;
nargmax = 8;
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

% Array position (m)
X0 = conf.X0;                    
Y0 = conf.Y0;
% Loudspeaker distcane
% NOTE: if dLS <= dx, than we have a continous loudspeaker
dLS = conf.LSdist;
%dLS = 0.00001;

% Reference position for the amplitude (correct reproduction of amplitude
% at y = yref).
yref = conf.yref;

% Use tapering window?
usetapwin = conf.usetapwin;

% xy resolution
xysamples = conf.xysamples;

% phase of omega
phase = conf.phase;

% Speed of sound
c = conf.c;


%% ===== Variables ======================================================

% General
omega = 2*pi*f;

% kx resolution (repetitions of the spectrum)
Nrep = 21;
kxsamples = Nrep * xysamples;
% Start and end index of the indexes in kx space that correspond to the
% given X vector
xstart = kxsamples/2 - (xysamples/2-1);
xend   = kxsamples/2 + xysamples/2;

% Geometry
y       = linspace(-0.1+Y0,Y+Y0,xysamples);
x       = linspace(-X/2,X/2,xysamples);
dx      = X/xysamples;  % x resolution
kxrep   = 2*pi / dx;
kx      = linspace(-kxrep/2,kxrep/2,kxsamples);

% indexes for evanescent contributions
idxpr = (( abs(kx) <= (omega/c) ));
idxev = (( abs(kx) >  (omega/c) ));


%% ===== Computation ====================================================

% Check if yref is in the given y space
if yref>max(y)
    error('%s: yref has be smaller than max(y) = %.2f',...
        upper(mfilename),max(y));
end


% Adjust the y position of the focused source due to the array position
%ys = ys + conf.Y0


%% ===== Spectrum of driving function =================================== 
%
%
% Dkx(x,omega) = 
%
%
Dkx = zeros(1,length(kx));

if strcmp('pw',src)
    
    % ===== PLANE WAVE ===================================================
    to_be_implemented(mfilename);

elseif strcmp('ps',src)
    
    % ===== POINT SOURCE =================================================
    to_be_implemented(mfilename);

elseif strcmp('fs',src)
    
    % ===== FOCUSED SOURCE ===============================================
    % D_2.5D using a point source (G_3D) as source model. The driving function
    % was calculated using the SDM method and then approximated by the large
    % argument approximation of the Hankel function (see Spors2009)
    %                ________                        ____________     
    %               |  y_ref |              / e^(i \|(w/c)^2-kx^2 ys), |kx|<|w/c|
    % Dkx(kx,w) = _ |--------  e^(i kx xs) <        ____________
    %              \|y_ref-ys               \ e^(-\|kx^2-(w/c)^2 ys),  |w/c|<|kx|
    %
    % NOTE: the phase term e^(i phase)  is only there in order to be able to
    % simulate different time steps
    %
    Dkx(idxpr) = exp(1i*kx(idxpr)*xs) .* ...
        exp(1i*sqrt( (omega/c).^2 - kx(idxpr).^2 )*ys) .* exp(1i*phase);
    Dkx(idxev) = exp(1i*kx(idxev)*xs) .* ...
        exp(-1*sqrt( kx(idxev).^2 - (omega/c).^2 )*ys) .* exp(1i*phase);

end


%% ===== Manipulation of the driving function ===========================
% The manipulations will be done in x space, therefore we have to first
% transform Dkx(kx,omega) to D(x,omega)

% Convert to x space
D = ifftx(Dkx,[],2);

% === Spatial aliasing ===
% Number of loudspeakers
nLS = number_of_loudspeaker(L,conf);
% Get the position of the loudspeakers
[LSpos,LSdir] = LSpos_linear(X0,Y0,(nLS-1)*dLS,nLS);
% Index the driving function D
idx = 1:length(D);
% Find the position of the first loudspeaker
% NOTE: we need regular index positions of the loudspeakers. Therefore we have
% to use a constant proportion from this position on. This means the lower your
% xysamples resolution the greater the error of the real loudspeaker placing.
idx1 = find(x>=LSpos(1,1),1,'first');
% Discretization by loudspeakers. Use only every dLS/dx point to simulate
% loudspeakers with a distance of dLS
didx = round(dLS/dx);
% Remove all indexes that contain a loudspeaker
idx = idx( mod(idx-idx1,didx)~=0 );
% Set the remaining values (all positions without a loudspeaker) to 0
D(idx) = 0;

% === Truncation ===
% Calculate the indices for the truncation from the relation of X and the
% array length L
% Relation of L and X
rel = L/X;
% Find center position of the array (conf.X0) in the x vector
% NOTE: D has to be a mirrored version of the x axis, therefore we are
% looking here for -X0
[a,X0idx] = find(x>=-conf.X0,1);
% Array position offset
Lstart = round(X0idx - xysamples/2 * rel);
Lend = round(X0idx + xysamples/2 * rel)-1;

% Check start and end values of the array and apply truncation
% Truncate the array in the left half of the kx space
% TODO: Check symetrie of Lsatrt and Lend!
if xstart+Lstart>2 && L~=Inf
    D(1:xstart+Lstart-2) = 0;
end
if xstart+Lend<=length(D) && L~=Inf
    D(xstart+Lend:end) = 0;
end

% === Tapering window ===
% Apply window
if(usetapwin)
    
    samples = Lend - (Lstart-1);
    onset = ceil(0.15*samples);
    offset = onset;
    win = hanningwin(onset,offset,samples);
    D(xstart+Lstart-1:xstart+Lend-1) = ...
        D(xstart+Lstart-1:xstart+Lend-1).*win';
    
end


% Transform back in kx space
Dkx = fftx(D,[],2);


%% ===== Spectrum of secondary sources ================================== 
Gkx = zeros(length(kx),length(y));

if strcmp(sourcetype,'ps')
    % point source
    [K,Y] = meshgrid(kx(idxpr),abs(y-Y0));
    Gkx(idxpr,:) = -1j/4 .* besselh(0,2,sqrt( (omega/c)^2 - K.^2 ).* Y)';
    [K,Y] = meshgrid(kx(idxev),abs(y-Y0));
    Gkx(idxev,:) = 1/(2*pi) .* besselk(0,sqrt( K.^2 - (omega/c)^2).* Y)';
elseif strcmp(sourcetype,'pw')
    % plane wave
    % TODO: implement me
    error('%s: The plane wave source type is not implemented at the moment', ...
        upper(mfilename));
else
    error(['%s: %s is a unknown source type. See help for avaiable ', ...
        'source types'],upper(mfilename));
end


%% ===== Reproduced wave field ========================================== 

% Pkx = Dkx * Gkx
Pkx = repmat(Dkx',1,length(y)) .* Gkx;


%% ===== Inverse spatial Fourier transformation ========================= 

P = ifftx( Pkx, [], 1);
P = P(xstart:xend,:);

% === Scale signal (at xs,yref) ===
% Find index
[a,xidx] = find(x>xs,1);
[a,yidx] = find(y>yref,1);
% Scale signal to 1
P = 1*P/abs(P(xidx,yidx));
