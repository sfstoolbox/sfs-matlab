function D = driving_function_mono_wfs_2d(x0,xs,src,f,conf)
%DRIVING_FUNCTION_MONO_WFS_2D returns the driving signal D for 2D WFS
%
%   Usage: D = driving_function_mono_wfs_2d(x0,xs,src,f,[conf])
%
%   Input parameters:
%       x0          - position and direction of the secondary source (m)
%       xs          - position of virtual source or direction of plane wave (m)
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs, ys are the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'ls' - line source
%                         'fs' - focused line source
%       f           - frequency of the monochromatic source (Hz)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal (1x1)
%
%   DRIVING_FUNCTION_MONO_WFS_2D(x0,xs,src,f,conf) returns the driving signal
%   for the given secondary source and desired source type (src).
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

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
isargsecondarysource(x0);
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
xref = position_vector(conf.xref);


%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain
%
% Omega
omega = 2*pi*f;

% Driving function D(x0,omega)

% Direction and position of secondary sources
nx0 = x0(4:6);
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
    % D_2D(x0,w) = 2i w/c n(xs) n(x0)  e^(-i w/c n(xs) x0)
    %
    % NOTE: the phase term e^(-i phase) is only there in order to be able to
    %       simulate different time steps
    %
    D = 2*1i*omega/c*nxs*nx0' * exp(-1i*omega/c*(nxs*x0')) * ...
        exp(-1i*phase);
    %

elseif strcmp('ps',src)

    % ===== POINT SOURCE =============================================
    %
    % D_2D using a point source
    %
    %                1  / iw      1    \  (x0-xs)nk
    % D_2D(x0,w) =  --- | -- - ------- |  --------- e^(-i w/c |x0-xs|)
    %               2pi \  c   |x0-xs| /  |x0-xs|^2
    %
    % NOTE: the phase term e^(-i phase) is only there in order to be able to
    %       simulate different time steps
    %
    D = 1/(2*pi) * ( (1i*omega)/c - 1/norm(x0-xs) ) * ...
        (x0-xs)*nx0' / norm(x0-xs)^2 * exp(-1i*omega/c*norm(x0-xs)) * ...
        exp(-1i*phase);
    %

elseif strcmp('ls',src)

    % ===== LINE SOURCE ==============================================
    %
    % D_2D using a line source
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

elseif strcmp('fs',src)

    % ===== FOCUSED SOURCE ===========================================
    %
    % D_2D using a line sink
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
else
    % No such source type for the driving function
    error('%s: src has to be one of "pw", "ps", "fs"!',upper(mfilename));
end
