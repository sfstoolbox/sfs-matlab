function [D] = driving_function_wfs_25d(x0,y0,phi,xs,ys,f,src,conf)
%DRIVING_FUNCTION_WFS_25D returns the driving signal D for 2.5D WFS
%   Usage: D = driving_function_wfs_25d(x0,y0,phi,xs,ys,f,src,conf)
%          D = driving_function_wfs_25d(x0,y0,phi,xs,ys,f,src)
%
%   Input parameters:
%       x0,y0,phi   - position and direction of the secondary source (m)
%       xs,ys       - position of point source or direction of plane wave (m)
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs, ys are the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal (1x1)
%
%   DRIVING_FUNCTION_WFS_25D(x0,y0,phi,xs,ys,f,src,conf) returns the driving signal
%   for the given secondary source and desired sourcei type (src).
%   The driving signal is calculated for the WFS 2.5 dimensional case in the
%   temporal domain.
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
nargmin = 8;
nargmax = 9;
error(nargchk(nargmin,nargmax,nargin));
isargscalar({x0,y0,phi,xs,ys},{'x0','y0','phi','xs','ys'});
isargpositivescalar({f},{'f'}),
isargchar({src},{'src'});
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct({conf},{'conf'});
end


%% ===== Configuration ==================================================

% phase of omega
phase = conf.phase;

% yref
yref = conf.yref;

% Plotting result
useplot = conf.useplot;

% Speed of sound
c = conf.c;


%% ===== Computation ====================================================

% Calculate the driving function in time-frequency domain
%
% Omega
omega = 2*pi*f;


% Driving function D(x0,omega)
% Activity of secondary sources
ls_activity = secondary_source_selection(x0,y0,phi,xs,ys,src);
if(ls_activity)
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
        D = -g0/(2*pi) * ( sqrt(omega/(1i*c)) + sqrt(1i*c/omega) / ...
            sqrt((x0-xs)^2+(y0-ys)^2) ) * ...
            (y0-ys)/((x0-xs)^2+(y0-ys)^2) * ...
            exp(-1i*omega/c*sqrt((x0-xs)^2+(y0-ys)^2)) * exp(-1i*phase);
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
        %D = -g0 * 1i*omega/(2*c) * ...
        %    (y0-ys)/(sqrt((x0-xs)^2+(y0-ys)^2) )^(3/2) * ...
        %    besselh(1,2,omega/c*sqrt((x0-xs)^2+(y0-ys)^2)) * exp(-1i*phase);
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
else
    D = 0;
end

% ===== Plotting =========================================================
if(useplot)
    figure; plot(D)
end
