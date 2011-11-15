function [D] = driving_function_mono_wfs_2d(x0,y0,phi,xs,ys,f,src,conf)
%DRIVING_FUNCTION_MONO_WFS_2D returns the driving signal D for 2D WFS
%   Usage: D = driving_function_mono_wfs_2d(x0,y0,phi,xs,ys,f,src,conf)
%          D = driving_function_mono_wfs_2d(x0,y0,phi,xs,ys,f,src)
%
%   Input parameters:
%       x0,y0,phi   - position and direction of the secondary source (m)
%       xs,ys       - position of virtual source or direction of plane wave (m)
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
%   DRIVING_FUNCTION_MONO_WFS_2D(x0,y0,phi,xs,ys,f,src,conf) returns the
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
nargmin = 7;
nargmax = 8;
error(nargchk(nargmin,nargmax,nargin));
isargscalar(x0,y0,phi,xs,ys);
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
% [xref yref]
xref = conf.xref;
yref = conf.yref;
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

    % Direction of secondary sources
    [nx0,ny0] = sph2cart(phi,0,1);


    if strcmp('pw',src)

        % ===== PLANE WAVE ===============================================
        % Use the position of the source as the direction vector for a plane
        % wave
        nxs = xs / sqrt(xs^2+ys^2);
        nys = ys / sqrt(xs^2+ys^2);
        %
        % ----------------------------------------------------------------
        % D_2D using a plane wave as source model
        %
        % D_2D(x0,w) = 2i w/c n(xs) n(x0)  e^(i w/c n(xs) x0)
        %
        % NOTE: the phase term e^(-i phase) is only there in order to be able to
        %       simulate different time steps
        %
        D = 2*1i*omega/c*[nxs nys]*[nx0 ny0]' * ...
            exp(1i*omega/c*(nxs*x0+nys*y0)) * ...
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
        D = 1/(2*pi) * ( (1i*omega)/c - 1/sqrt((x0-xs)^2+(y0-ys)^2) ) * ...
            ([x0 y0]-[xs ys])*[nx0 ny0]' / ((x0-xs)^2+(y0-ys)^2) * ...
            exp(1i*omega/c*sqrt((x0-xs)^2+(y0-ys)^2)) * exp(-1i*phase);
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
        D = -1i*omega/(2*c) * ...
            ([x0 y0]-[xs ys])*[nx0 ny0]' / ...
            (sqrt((x0-xs)^2+(y0-ys)^2) )^(3/2) * ...
            besselh(1,1,omega/c*sqrt((x0-xs)^2+(y0-ys)^2)) * exp(-1i*phase);
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
        D = -1i*omega/(2*c) * ...
            ([x0 y0]-[xs ys])*[nx0 ny0]' / ...
            (sqrt((x0-xs)^2+(y0-ys)^2) )^(3/2) * ...
            besselh(1,2,omega/c*sqrt((x0-xs)^2+(y0-ys)^2)) * exp(-1i*phase);
        %
    else
        % No such source type for the driving function
        error('%s: src has to be one of "pw", "ps", "fs"!',upper(mfilename));
    end
else
    D = 0;
end
