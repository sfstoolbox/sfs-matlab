function [f,S] = freq_WFS_25D(X,Y,xs,ys,L,src,conf)
%FREQ_WFS_25D simulates the frequency response of a wave field for 2.5D WFS 
%   Usage: [f,S] = freq_WFS_25D(X,Y,xs,ys,L,src,conf)
%          [f,S] = freq_WFS_25D(X,Y,xs,ys,L,src)
%
%   Input parameters:
%       X           - x-axis position of the frequency response (m)
%       Y           - y-axis position of the frequency response (m)
%       xs          - x position of point source (m)
%       ys          - y position of point source (m)
%       L           - array length (m)
%       src         - source type of the virtual source
%                         'pw' -plane wave
%                         'ps' - point source
%                         'fs' - focused source
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       f           - corresponding frequency (x) axis
%       S           - simulated frequency response
%
%   FREQ_WFS_25D(X,Y,xs,ys,L,src,conf) simulates the frequency response of
%   the wave field at the given position [X,Y]. The wave field is simulated for 
%   the given source type (src) using a WFS 2.5 dimensional driving function in
%   the temporal domain. 
%
%   References:
%       Spors2009 - Physical and Perceptual Properties of Focused Sources in
%           Wave Field Synthesis (AES127)
%       Spors2010 - Analysis and Improvement of Pre-equalization in
%           2.5-Dimensional Wave Field Synthesis (AES128)
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: wf_WFS_25D

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 6;
nargmax = 7;
error(nargchk(nargmin,nargmax,nargin));

if ~isnumeric(X) || ~isscalar(X)
    error('%s: X has to be a scalar!',upper(mfilename));
end
if ~isnumeric(Y) || ~isscalar(Y)
    error('%s: Y has to be a scalar!',upper(mfilename));
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
    conf = SFS_config;
end

% Array position (m)
X0 = conf.X0;                    
Y0 = conf.Y0;

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

% Phase of the wave field
phase = conf.phase;

% Plotting result
useplot = conf.useplot;

% Speed of sound
c = conf.c;


%% ===== Computation ====================================================
% Calculate the wave field in time domain 
%
% Number of loudspeakers
nLS = number_of_loudspeaker(L,conf);
% Get the position of the loudspeakers
[LSpos,LSdir] = LSpos_linear(X0,Y0,(nLS-1)*dLS,nLS);
% The array loudspeakers are all positioned at Y0 in the y direction 
% (linear array)
y0 = Y0;
% Generate frequencies (10^0-10^5)
f = logspace(0,5,500);
% We want only frequencies until f = 20000Hz
idx = find(f>20000,1);
f = f(1:idx);
S = zeros(1,length(f));
% Tapering window
win = tapwin(L,conf);
% Get the result for all frequencies
for ii = 1:length(f)
    omega = 2*pi*f(ii);
    P = 0;
    % Iteration index for tapwin
    jj = 1;
    % Integration over secondary source positions
    for x0 = LSpos(1,1):dLS:LSpos(1,end)
        
        % ================================================================
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
        
        % ================================================================
        % Driving function D(x0,omega)
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
        % NOTE: this driving function comes from an error in the calculation
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
        %warning('%s: use of a non standard driving function!',upper(mfilename);
        %D = -g0/(2*pi) * ( sqrt(omega/(1i*c)) + sqrt(1i*c/omega) ) * ...
        %    (y0-ys)/(sqrt((x0-xs)^2+(y0-ys)^2) )^(3/2) * ...
        %    exp(-1i*omega/c*sqrt((x0-xs)^2+(y0-ys)^2)) * exp(-1i*phase);
        
        % --------------------------------------------------------------------
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
        %warning('%s: use of a non standard driving function!',upper(mfilename);
        %D = -g0 * 1i*omega/(2*c) * ...
        %    (y0-ys)/(sqrt((x0-xs)^2+(y0-ys)^2) )^(3/2) * ...
        %    besselh(1,2,omega/c*sqrt((x0-xs)^2+(y0-ys)^2)) * ...
        %    exp(-1i*phase);
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
        %warning('%s: use of a non standard driving function!',upper(mfilename));
        D = -g0/(2*pi) * sqrt(1i*omega/c) * ...
            (y0-ys)/(sqrt((x0-xs)^2+(y0-ys)^2))^(3/2) * ...
            exp(-1i*omega/c*sqrt((x0-xs)^2+(y0-ys)^2)) * exp(-1i*phase);
        %
        % --------------------------------------------------------------------
        % Integration
        %              /
        % P(x,omega) = | D(x0,omega) G(x-xs,omega) dx0
        %              /
        % 
        % see: Spors2009, Williams1993 p. 36
        %
        P = P + win(jj)*D.*G;
        jj = jj+1;
    end
    S(ii) = abs(P);
end


% ===== Plotting =========================================================
if(useplot)
    figure; semilogx(f,20*log10(S));
    ylabel('Amplitude (dB)');
    xlabel('Frequency (Hz)');
end
