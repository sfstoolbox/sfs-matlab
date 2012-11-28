function [D] = driving_function_mono_wfs_25d(x0,xs,src,f,conf)
%DRIVING_FUNCTION_MONO_WFS_25D returns the driving signal D for 2.5D WFS
%
%   Usage: D = driving_function_mono_wfs_25d(x0,xs,src,f,[conf])
%
%   Input parameters:
%       x0          - position and direction of the secondary source (m)
%       xs          - position of virtual source or direction of plane wave (m)
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs is the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       f           - frequency of the monochromatic source (Hz)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal (1x1)
%
%   DRIVING_FUNCTION_MONO_WFS_25D(x0,xs,f,src,conf) returns the
%   driving signal for the given secondary source and desired source type (src).
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
% xref
xref = position_vector(conf.xref);
% Speed of sound
c = conf.c;


%% ===== Computation ====================================================

% Calculate the driving function in time-frequency domain
%
% Omega
omega = 2*pi*f;

% Driving function D(x0,omega)

% Direction and position of secondary sources
nx0 = x0(:,4:6);
x0 = x0(:,1:3);

% Constant amplitude factor g0
g0 = sqrt(2*pi*norm(xref-x0));

% Initialize empty driving function
D = zeros(size(x0,1));

if strcmp('pw',src)

    % ===== PLANE WAVE ===============================================
    % Use the position of the source as the direction vector for a plane
    % wave
    nxs = xs / norm(xs);
    %
    % ----------------------------------------------------------------
    % D_25D using a plane wave as source model
    %                                   ___
    %                                  | w |
    % D_25D(x0,w) = 2 g0 n(xs) n(x0) _ |---  e^(-i w/c n(xs) x0)
    %                                 \|i c
    %
    % NOTE: the phase term e^(-i phase) is only there in order to be able to
    %       simulate different time steps
    %
    D = 2*g0* nxs*nx0' * sqrt(omega/(1i*c)) * ...
        exp(-1i*omega/c*(nxs*x0')) * ...
        exp(-1i*phase);


elseif strcmp('ps',src)

    % ===== POINT SOURCE =============================================
    %
    % ----------------------------------------------------------------
    % D_25D using a point source as source model
    %
    % D_25D(x0,w) =
    %             ___       ___
    %    g0  /   | w |     |i c|    1    \  (x0-xs)nk
    %   ---  | _ |---  - _ |---  ------- |  --------- e^(-i w/c |x0-xs|)
    %   2pi  \  \|i c     \| w   |x0-xs| /  |x0-xs|^2
    %
    % NOTE: the phase term e^(-i phase) is only there in order to be able to
    %       simulate different time steps
    %
    D = g0/(2*pi) * ( sqrt(omega/(1i*c)) - sqrt(1i*c/omega) / ...
        norm(x0-xs) ) * (x0-xs)*nx0' / norm(x0-xs)^2 * ...
        exp(-1i*omega/c*norm(x0-xs)) * exp(-1i*phase);

elseif strcmp('fs',src)

    % ===== FOCUSED SOURCE ===========================================
    %
    % ----------------------------------------------------------------
    % D_25D using a point sink as source model
    %
    % D_25D(x0,w) =
    %             ___       ___
    %   -g0  /   | w |     |i c|    1    \  (x0-xs)nk
    %   ---  | _ |---  + _ |---  ------- |  --------- e^(i w/c |x0-xs|)
    %   2pi  \  \|i c     \| w   |x0-xs| /  |x0-xs|^2
    %
    % NOTE: the phase term e^(-i phase) is only there in order to be able to
    %       simulate different time steps
    %
    %D = -g0/(2*pi) * ( sqrt(omega/(1i*c)) + sqrt(1i*c/omega) / ...
    %    norm(x0-xs) ) * (x0-xs)*nx0' / norm(x0-xs)^2 * ...
    %    exp(1i*omega/c*norm(x0-xs)) * exp(-1i*phase);
    %
    % ----------------------------------------------------------------
    % Alternative Driving Functions for a focused source:
    %
    % D_25D using a line sink with point source amplitude characteristics as
    % source (see Spors2009).
    %
    %                   iw (x0-xs)nk   (1)/ w         \
    % D_25D(x0,w) = -g0 -- --------- H1  | - |x0-xs| |
    %                   2c |x0-xs|        \ c         /
    %
    % NOTE: the phase term e^(-i phase) is only there in order to be able to
    %       simulate different time steps
    %
    %D = -g0 * 1i*omega/(2*c) * (x0-xs)*nx0' / norm(x0-xs)^(3/2) * ...
    %    besselh(1,1,omega/c*norm(x0-xs)) * exp(-1i*phase);
    %
    % --------------------------------------------------------------------
    % D_25D using a line sink with point amplitue characteristic as source
    % and the large argument approximation of the driving function above.
    % This results in the "traditional" driving function, derived in
    % Verheijen1997 (see Spors2009).
    %                     _____
    %                    |i w  |   (x0-xs)nk
    % D_25D(x0,w) = g0 _ |-----  ------------- e^(i w/c |x0-xs|)
    %                   \|2pi c  |x0-xs|^(3/2)
    %
    % NOTE: the phase term e^(-i phase) is only there in order to be able to
    %       simulate different time steps
    %
    D = g0 * sqrt(1i*omega/(2*pi*c)) * ...
        (x0-xs)*nx0' / norm(x0-xs)^(3/2) * ...
        exp(1i*omega/c*norm(x0-xs)) * exp(-1i*phase);
    %
else
    % No such source type for the driving function
    error('%s: src has to be one of "pw", "ps", "fs"!',upper(mfilename));
end
