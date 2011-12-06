function [D] = driving_function_mono_wfs_2d(x0,xs,f,src,conf)
%DRIVING_FUNCTION_MONO_WFS_2D returns the driving signal D for 2D WFS
%   Usage: D = driving_function_mono_wfs_2d(x0,xs,f,src,conf)
%          D = driving_function_mono_wfs_2d(x0,xs,f,src)
%
%   Input parameters:
%       x0          - position and direction of the secondary source (m)
%       xs          - position of virtual source or direction of plane wave (m)
%       f           - frequency of the monochromatic source (Hz)
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs, ys are the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'ls' - line source
%                         'fs' - focused line source
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal (1x1)
%
%   DRIVING_FUNCTION_MONO_WFS_2D(x0,xs,f,src,conf) returns the
%   driving signal for the given secondary source and desired source type (src).
%   The driving signal is calculated for the WFS 2 dimensional case in the
%   temporal domain.
%
%   References:
%       Spors2009 - Physical and Perceptual Properties of Focused Sources in
%           Wave Field Synthesis (AES127)
%       Spors2010 - Analysis and Improvement of Pre-equalization in
%           2.5-Dimensional Wave Field Synthesis (AES128)
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: plot_wavefield, wave_field_mono_wfs_25d,
%             driving_function_imp_wfs_25d

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 5;
error(nargchk(nargmin,nargmax,nargin));
isargsecondarysource(x0);
isargposition(xs);
xs = position_vector(xs);
isargpositivescalar(f);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================

% phase of omega
phase = conf.phase;
% Speed of sound
c = conf.c;


%% ===== Computation ====================================================

% Calculate the driving function in time-frequency domain
%
% Omega
omega = 2*pi*f;


% Driving function D(x0,omega)
% Activity of secondary sources
ls_activity = secondary_source_selection(x0,xs,src);
if(ls_activity)

    % Direction and position of secondary sources
    nx0 = secondary_source_direction(x0);
    x0 = x0(1:3);


    if strcmp('pw',src)

        % ===== PLANE WAVE ===============================================
        % Use the position of the source as the direction vector for a plane
        % wave
        nxs = xs / norm(xs);
        %
        % ----------------------------------------------------------------
        % D_2D using a plane wave as source model
        %
        % D_2D(x0,w) = 2i w/c n(xs) n(x0)  e^(i w/c n(xs) x0)
        %
        % NOTE: the phase term e^(-i phase) is only there in order to be able to
        %       simulate different time steps
        %
        D = 2*1i*omega/c*nxs*nx0' * exp(1i*omega/c*(nxs*x0')) * ...
            exp(-1i*phase);
        %

    elseif strcmp('ps',src)

        % ===== POINT SOURCE =============================================
        %
        % D_2D using a point source
        %
        %                1  / iw      1    \  (x0-xs)nk
        % D_2D(x0,w) =  --- | -- - ------- |  --------- e^(i w/c |x0-xs|)
        %               2pi \  c   |x0-xs| /  |x0-xs|^2
        %
        % NOTE: the phase term e^(-i phase) is only there in order to be able to
        %       simulate different time steps
        %
        D = 1/(2*pi) * ( (1i*omega)/c - 1/norm(x0-xs) ) * ...
            (x0-xs)*nx0' / norm(x0-xs)^2 * exp(1i*omega/c*norm(x0-xs)) * ...
            exp(-1i*phase);
        %

    elseif strcmp('ls',src)

        % ===== LINE SOURCE ==============================================
        %
        % D_2D using a line source
        %
        %                 iw (x0-xs)nk  (1)/ w         \
        % D_2D(x0,w) =  - -- --------- H1  | - |x0-xs| |
        %                 2c  |x0-xs|      \ c         /
        %
        % NOTE: the phase term e^(-i phase) is only there in order to be able to
        %       simulate different time steps
        %
        D = -1i*omega/(2*c) * (x0-xs)*nx0' / norm(x0-xs)^(3/2) * ...
            besselh(1,1,omega/c*norm(x0-xs)) * exp(-1i*phase);
        %

    elseif strcmp('fs',src)

        % ===== FOCUSED SOURCE ===========================================
        %
        % D_2D using a line sink
        %
        %                 iw (x0-xs)nk  (2)/ w         \
        % D_2D(x0,w) =  - -- --------- H1  | - |x0-xs| |
        %                 2c  |x0-xs|      \ c         /
        %
        % NOTE: the phase term e^(-i phase) is only there in order to be able to
        %       simulate different time steps
        %
        D = -1i*omega/(2*c) * (x0-xs)*nx0' / norm(x0-xs)^(3/2) * ...
            besselh(1,2,omega/c*norm(x0-xs)) * exp(-1i*phase);
        %
    else
        % No such source type for the driving function
        error('%s: src has to be one of "pw", "ps", "fs"!',upper(mfilename));
    end
else
    D = 0;
end
