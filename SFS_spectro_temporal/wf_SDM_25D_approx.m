function [x,y,P] = wf_SDM_25D_approx(X,Y,xs,ys,L,f,src,conf)
%WF_SDM_25D_APPROX simulates the wave field of a given source for 2.5D SDM
%   Usage: [x,y,P] = wf_SDM_25D_approx(X,Y,xs,ys,L,f,src,conf)
%          [x,y,P] = wf_SDM_25D_approx(X,Y,xs,ys,L,f,src)
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
%       P           - Simulated wave field in the format P = P(x,y)
%
%   WF_SDM_25D_approx(X,Y,xs,ys,L,f,src,conf) simulates a wave field of the given
%   source type (src) using the large argument approximation of a SDM 2.5 
%   dimensional driving function in the spectro-temporal domain. 
%   To plot the result use plot_wavefield(x,y,P).
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

% Using large argument approximation (=> driving function as in the WFS case,
% but with an amplitude factor)
large_argument_approximation = 1;

% Array position (m)
X0 = conf.X0;                    
Y0 = conf.Y0;
% Loudspeaker distcane
% NOTE: if dLS <= dx, than we have a continous loudspeaker
dLS = conf.LSdist;

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
    if(large_argument_approximation)
        % Large argument approximation for Bessel function
        %                              ______
        %                             | yref | 
        % D_25D(kx,w) = e^(i kx xs) _ |-------  ...
        %                            \|yref-ys  
        %
        %                   / e^(i \|(w/c)^2-kx^2 ys), |kx|<|w/c|
        %                  <
        %                   \ e^(-\|kx^2-(w/c)^2 ys),  |kx|>|w/c|
        %
        % NOTE: the phase term e^(i phase) is only there in order to be able
        %       to simulate different time steps
        %
        Dkx(idxpr) = exp(1i*kx(idxpr)*xs) .* sqrt(yref/(yref-ys)) .* ...
            exp(1i*sqrt( (omega/c).^2 - kx(idxpr).^2 )*ys) .* ...
            exp(1i*phase);
        if(withev)
            Dkx(idxev) = exp(1i*kx(idxev)*xs) .* sqrt(yref/(yref-ys)) .* ...
                exp(-1*sqrt( kx(idxev).^2 - (omega/c).^2 )*ys) .* ...
                exp(1i*phase);
        end
    else
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
    end

else
    % No such source type for the driving function
    error('%s: src has to be one of "pw", "ps", "fs"!',upper(mfilename));
end


%% ===== Manipulation of the driving function ===========================
% The manipulations will be done in x space, therefore we have to first
% transform Dkx(kx,omega) to D(x,omega)

% Convert to x space
D = ifftx(Dkx,[],2);


% === Spatial aliasing ===
% Discretization by loudspeakers. Use only every dLS/dx point to simulate
% loudspeakers with a distance of dLS
idx = 1:length(D);
loudspeaker_idx = ceil(dLS/dx);
% Remove all indexes that contain a loudspeaker
idx = idx( mod(idx-1, loudspeaker_idx)~=0 );
% Set the remaining values to 0
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

P = ifftx( Pkx, [], 1);
P = P(xstart:xend,:);

% === Scale signal (at xs,yref) ===
% Find index
[a,xidx] = find(x>xs,1);
[a,yidx] = find(y>yref,1);
% Scale signal to 1
P = 1*P/abs(P(xidx,yidx));
