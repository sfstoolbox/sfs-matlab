function [x,y,P] = wave_field_mono_sdm_25d(X,Y,xs,src,f,L,conf)
%WAVE_FIELD_MONO_SDM_25D simulates a wave field for 2.5D NFC-HOA
%
%   Usage: [x,y,P] = wave_field_mono_sdm_25d(X,Y,xs,src,f,L,[conf])
%
%   Input parameters:
%       X           - [xmin,xmax]
%       Y           - [ymin,ymax]
%       xs          - position of point source (m)
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs is the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       f           - monochromatic frequency (Hz)
%       L           - array length (m)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       x           - corresponding x axis
%       y           - corresponding y axis
%       P           - Simulated wave field
%
%   WAVE_FIELD_MONO_SDM_25D(X,Y,xs,src,f,L,conf) simulates a wave
%   field of the given source type (src) using a SDM 2.5 dimensional driving
%   function in the space/time-frequency domain. This means by calculating 
%   the integral for P with a summation.
%   To plot the result use plot_wavefield(x,y,P).
%
%   References:
%       
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: plot_wavefield, wave_field_mono_sdm_25d_kx

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
nargmin = 6;
nargmax = 7;
narginchk(nargmin,nargmax);
isargvector(X,Y);
xs = position_vector(xs);
isargpositivescalar(L,f);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
% Plotting result
useplot = conf.useplot;


%% ===== Computation ====================================================
% Calculate the wave field in time-frequency domain
%
% Get the position of the loudspeakers
x0 = secondary_source_positions(L,conf);
% Create a x-y-grid to avoid a loop
[xx,yy,x,y] = xy_grid(X,Y,conf);
% Initialize empty wave field
P = zeros(length(y),length(x));
% Driving function D(x0,omega)
D = driving_function_mono_sdm_25d(x0,xs,src,f,conf);
% Generate tapering window
win = tapering_window(x0,conf);

% Integration over secondary source positions
for ii = 1:size(x0,1)

    % ====================================================================
    % Secondary source model G(x-x0,omega)
    % This is the model for the loudspeakers we apply. We use closed cabinet
    % loudspeakers and therefore point sources.
    G = point_source(xx,yy,x0(ii,1:3),f);

    % ====================================================================
    % Integration
    %              /
    % P(x,omega) = | D(x0,omega) G(x-x0,omega) dx0
    %              /
    %
    % see: Spors2009, Williams1993 p. 36
    P = P + win(ii).*D(ii).*G;

end

% === Scale signal (at xref) ===
P = norm_wave_field(P,x,y,conf);


% ===== Plotting =========================================================
if(useplot)
    plot_wavefield(x,y,P,x0,win,conf);
end
