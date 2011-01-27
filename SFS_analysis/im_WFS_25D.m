function [x,y,P] = im_WFS_25D(X,Y,xs,ys,L,f,src,conf)
%WF_WFS_25D simulates the wave field of a given source for 2.5D WFS 
%   Usage: [x,y,P] = wf_WFS_25D(X,Y,xs,ys,L,f,src,conf)
%          [x,y,P] = wf_WFS_25D(X,Y,xs,ys,L,f,src)
%
%   Input parameters:
%       X           - length of the X axis (m) [xaxis: -X/2:X/2]
%       Y           - length of the Y axis (m) [yaxis: -0.1:Y]
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
%   WF_WFS_25D(X,Y,xs,ys,L,f,src,conf) simulates a wave field of the given 
%   source type (src) using a WFS 2.5 dimensional driving function in the 
%   temporal domain. This means by calculating the integral for P with a 
%   summation. 
%   To plot the result use plot_wavefield(x,y,P).
%
%   References:
%       Spors2009 - Physical and Perceptual Properties of Focused Sources in
%           Wave Field Synthesis (AES127)
%       Spors2010 - Analysis and Improvement of Pre-equalization in
%           2.5-Dimensional Wave Field Synthesis (AES128)
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: plot_wavefield, wf_SDM_25D

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
% Check if the focused source is positioned before the loudspeaker array.
if strcmp('fs',src) && ys<=Y0
    error('%s: ys has to be greater than Y0 for a focused source.', ...
        upper(mfilename));
elseif strcmp('ps',src) && ys>=Y0
    error('%s: ys has to be smaller than Y0 for a point source.', ...
        upper(mfilename));
end

% Loudspeaker distcane
% NOTE: if dLS <= dx, than we have a continous loudspeaker
dLS = conf.LSdist;

% Reference position for the amplitude (correct reproduction of amplitude
% at y = yref).
yref = conf.yref;

% Use tapering window?
usetapwin = conf.usetapwin;

% xy resolution
xysamples = conf.xysamples;

% phase of omega
phase = conf.phase;

% Plotting result
useplot = conf.useplot;

% Speed of sound
c = conf.c;

% Sampling rate
fs = conf.fs;


%% ===== Variables ======================================================

% General
omega = 2*pi*f;

% Geometry
y = linspace(-0.1+Y0,Y+Y0,xysamples);
x = linspace(-X/2,X/2,xysamples);


%% ===== Computation ====================================================

% Check if yref is in the given y space
if yref>max(y)
    error('%s: yref has be smaller than max(y) = %.2f',...
        upper(mfilename),max(y));
end


% Calculate the wave field in time-frequency domain 
%
% Number of loudspeakers
nLS = number_of_loudspeaker(L,conf);
% Get the position of the loudspeakers
[LSpos,LSdir] = LSpos_linear(X0,Y0,(nLS-1)*dLS,nLS);
% Create a x-y-grid to avoid a loop
[Y,X] = meshgrid(y,x);
f = 1:10;
% Initialize empty wave field
P = zeros(length(x),length(y),length(f));
% The array loudspeakers are all positioned at Y0 in the y direction 
% (linear array)
y0 = Y0;
% Generate tapering window 
% NOTE: if you have disabled tapering window, this will give you back ones()
win = tapwin(L,conf);
% Iteration over all frequencies
for ff = 1:length(f)
    
    omega = 2*pi*f(ff);

    % Iteration index for tapwin
    ii = 1;
    
    % Integration over secondary source positions
    for x0 = LSpos(1,1):dLS:LSpos(1,end)

        % ====================================================================
        % Secondary source model
        % This is the model for the loudspeakers we apply. We use closed cabinet
        % loudspeakers and therefore the 3D Green's function is our model.
        %
        %              1  e^(-i w/c |x-xs|)
        % G(x-xs,w) = --- -----------------
        %             4pi      |x-xs|
        %
        % see: Williams1999, p. 198
        %
        G = 1/(4*pi) * exp(i*omega/c.*sqrt((X-x0).^2+(Y-y0).^2)) ./ ...
                sqrt((X-x0).^2+(Y-y0).^2);
        
        % ====================================================================
        % Driving function D(x0,omega)
        %
        if strcmp('pw',src)

            % ===== PLANE WAVE ===============================================
            % Constant pre-equalization factor g0
            g0 = sqrt(2*pi*abs(yref-y0));
            % Use the position of the source as the direction vector for a plane
            % wave
            nx0 = xs / sqrt(xs^2+ys^2);
            ny0 = ys / sqrt(xs^2+ys^2);
            %
            % ----------------------------------------------------------------
            % D_25D using a plane wave as source model
            %                             ___
            %                            | w |
            % D_25D(x0,w) = 2 g0 n_ky0 _ |---  e^(i w/c nk x0)  
            %                           \|i c
            %
            % NOTE: the phase term e^(-i phase) is only there in order to be able to
            %       simulate different time steps
            %
            D = 2*g0*ny0*sqrt(omega/(1i*c)) * e^(1i*omega/c*(nx0*x0+ny0*y0)) * ...
                exp(-1i*phase);
            

        elseif strcmp('ps',src)

            % ===== POINT SOURCE =============================================
            % Constant pre-equalization factor g0
            g0 = sqrt(2*pi*abs(yref-y0));
            %
            % ----------------------------------------------------------------
            % D_25D using a point source as source model
            %
            % D_25D(x0,w) = 
            %             ___       ___
            %    g0  /   | w |     |i c|    1    \    y0-ys  
            %   ---  | _ |---  - _ |---  ------- |  --------- e^(i w/c |x0-xs|)
            %   2pi  \  \|i c     \| w   |x0-xs| /  |x0-xs|^2
            %
            % NOTE: the phase term e^(-i phase) is only there in order to be able to
            %       simulate different time steps
            %
            D = g0/(2*pi) * ( sqrt(omega/(1i*c)) - sqrt(1i*c/omega) / ...
                sqrt((x0-xs)^2+(y0-ys)^2) ) * ...
                (y0-ys)/((x0-xs)^2+(y0-ys)^2) * ...
                exp(1i*omega/c*sqrt((x0-xs)^2+(y0-ys)^2)) * exp(-1i*phase);
            % ----------------------------------------------------------------
            % D_25D using a point source as source (G_3D)
            % NOTE: hier habe ich mich verrechnet!
            %
            % D_25D(x0,w) = 
            %              ___       ___
            %    -g0  /   | w |     |i c| \      y0-ys
            %    ---  | _ |---  + _ |---  |  ------------  e^(i w/c |x0-xs|)
            %    2pi  \  \|i c     \| w   /  |x0-xs|^(3/2)
            %    
            % NOTE: the phase term e^(-i phase) is only there in order to be able to
            %       simulate different time steps
            %
            %D = -g0/(2*pi) * ( sqrt(omega/(1i*c)) + sqrt(1i*c/omega) ) * ...
            %    (y0-ys)/(sqrt((x0-xs)^2+(y0-ys)^2) )^(3/2) * ...
            %    exp(1i*omega/c*sqrt((x0-xs)^2+(y0-ys)^2)) * exp(-1i*phase);
            

        elseif strcmp('fs',src)

            % ===== FOCUSED SOURCE ===========================================
            % 
            % Constant pre-equalization factor g0 (see Spors2009)
            g0 = sqrt(2*pi*abs(yref-y0));
            %
            % ----------------------------------------------------------------
            % D_25D using a point sink as source model
            %
            % D_25D(x0,w) = 
            %             ___       ___
            %   -g0  /   | w |     |i c|    1    \    y0-ys
            %   ---  | _ |---  + _ |---  ------- |  --------- e^(-i w/c |x0-xs|)
            %   2pi  \  \|i c     \| w   |x0-xs| /  |x0-xs|^2
            %
            % NOTE: the phase term e^(-i phase) is only there in order to be able to
            %       simulate different time steps
            %
            %D = -g0/(2*pi) * ( sqrt(omega/(1i*c)) + sqrt(1i*c/omega) / ...
            %    sqrt((x0-xs)^2+(y0-ys)^2) ) * ...
            %    (y0-ys)/((x0-xs)^2+(y0-ys)^2) * ...
            %    exp(-1i*omega/c*sqrt((x0-xs)^2+(y0-ys)^2)) * exp(-1i*phase);
            %
            % ----------------------------------------------------------------
            % D_25D using a point sink as source (G_3D)
            % NOTE: hier habe ich mich verrechnet
            %
            % D_25D(x0,w) = 
            %              ___       ___
            %    -g0  /   | w |     |i c| \      y0-ys
            %    ---  | _ |---  + _ |---  |  ------------  e^(-i w/c |x0-xs|)
            %    2pi  \  \|i c     \| w   /  |x0-xs|^(3/2)
            %    
            % NOTE: the phase term e^(-i phase) is only there in order to be able to
            %       simulate different time steps
            %
            %D = -g0/(2*pi) * ( sqrt(omega/(1i*c)) + sqrt(1i*c/omega) ) * ...
            %    (y0-ys)/(sqrt((x0-xs)^2+(y0-ys)^2) )^(3/2) * ...
            %    exp(-1i*omega/c*sqrt((x0-xs)^2+(y0-ys)^2)) * exp(-1i*phase);
            %
            % ----------------------------------------------------------------
            % Alternative Driving Functions for a focused source:
            %
            % D_25D using a line sink with point source amplitude characteristics as
            % source (see Spors2009).
            %
            %                   iw  y0-ys    (2)/ w         \
            % D_25D(x0,w) = -g0 -- -------  H1  | - |x0-xs| |  
            %                   2c |x0-xs|      \ c         /
            %
            % NOTE: the phase term e^(-i phase) is only there in order to be able to
            %       simulate different time steps
            %
            D = -g0 * 1i*omega/(2*c) * ...
                (y0-ys)/(sqrt((x0-xs)^2+(y0-ys)^2) )^(3/2) * ...
                besselh(1,2,omega/c*sqrt((x0-xs)^2+(y0-ys)^2)) * exp(-1i*phase);
            %
            % --------------------------------------------------------------------
            % D_25D using a line sink with point amplitue characteristic as source
            % and the large argument approximation of the driving function above. 
            % This results in the "traditional" driving function, derived in
            % Verheijen1997 (see Spors2009).
            %                        ___
            %                 g0    |i w|    y0-ys
            % D_25D(x0,w) = - --- _ |---  ------------- e^(-i w/c |x0-xs|)
            %                 2pi  \| c   |x0-xs|^(3/2)
            %
            % NOTE: the phase term e^(-i phase) is only there in order to be able to
            %       simulate different time steps
            %
            %D = -g0/(2*pi) * sqrt(1i*omega/c) * ...
            %    (y0-ys)/(sqrt((x0-xs)^2+(y0-ys)^2))^(3/2) * ...
            %    exp(-1i*omega/c*sqrt((x0-xs)^2+(y0-ys)^2)) * exp(-1i*phase);
            %
        else
            % No such source type for the driving function
            error('%s: src has to be one of "pw", "ps", "fs"!',upper(mfilename));
        end
        
        % --------------------------------------------------------------------
        % Integration
        %              /
        % P(x,omega) = | D(x0,omega) G(x-xs,omega) dx0
        %              /
        % 
        % see: Spors2009, Williams1993 p. 36
        % 
        % NOTE: win(ii) is the factor of the tapering window in order to have fewer
        % truncation artifacts. If you don't use a tapering window win(ii) will
        % always be one.
        P(:,:,ff) = P(:,:,ff) + win(ii)*D.*G;

        ii = ii+1;

    end

end

% Perform a IFFT
% NOTE: this should result in a real signal. If not something has gone wrong!
% The first entry of the spectrum has no phase and is just a sum over the time
% signal
am1 = length(f) * ones(length(x),length(y),length(f)+1);
ph1 = zeros(length(x),length(y),length(f)+1);
am1(:,:,2:end) = abs(P);
ph1(:,:,2:end) = angle(P);
am2 = am1(:,:,end-1:-1:2);
ph2 = -1*ph1(:,:,end-1:-1:2);
am = cat(3,am1,am2);
ph = cat(3,ph1,ph2);
spec = am .* exp(1i*ph);
p = ifft(spec,[],3);


% === Scale signal (at xs,yref) ===
% Find index
%[a,xidx] = find(x>xs,1);
%[a,yidx] = find(y>yref,1);
% Scale signal to 1
%P = 1*P/abs(P(xidx,yidx));


