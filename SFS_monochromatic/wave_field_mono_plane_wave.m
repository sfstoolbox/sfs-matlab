function [x,y,P] = wave_field_mono_plane_wave(X,Y,xs,f,conf)
%WAVE_FIELD_MONO_PLANE_WAVE simulates a wave field of a plane wave
%
%   Usage: [x,y,P] = wave_field_mono_plane_wave(X,Y,xs,f,[conf])
%
%   Input parameters:
%       X           - [xmin,xmax]
%       Y           - [ymin,ymax]
%       xs          - direction of the plane wave
%       f           - monochromatic frequency (Hz)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       x           - corresponding x axis
%       y           - corresponding y axis
%       P           - Simulated wave field
%
%   WAVE_FIELD_MONO_PLANE_WAVE(X,Y,xs,f,conf) simulates a wave
%   field of a plane wave going in the direction xs.
%   To plot the result use plot_wavefield(x,y,P).
%
%   References:
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: plot_wavefield, wave_field_mono_point_source

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
isargvector(X,Y);
xs = position_vector(xs);
isargpositivescalar(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
% Plotting result
useplot = conf.useplot;


%% ===== Computation ====================================================
% Create a x-y-grid to avoid a loop
[xx,yy,x,y] = xy_grid(X,Y,conf);
% Source model for a point source G(x,omega)
P = plane_wave(xx,yy,xs,f);
% Scale signal (at xref)
P = norm_wave_field(P,x,y,conf);


% ===== Plotting =========================================================
if(useplot)
    plot_wavefield(x,y,P,conf);
end
