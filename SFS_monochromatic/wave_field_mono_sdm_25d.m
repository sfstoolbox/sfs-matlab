function [x,y,z,P,x0] = wave_field_mono_sdm_25d(X,Y,Z,xs,src,f,L,conf)
%WAVE_FIELD_MONO_SDM_25D simulates a wave field for 2.5D NFC-HOA
%
%   Usage: [x,y,P] = wave_field_mono_sdm_25d(X,Y,Z,xs,src,f,L,[conf])
%
%   Input parameters:
%       X           - [xmin,xmax]
%       Y           - [ymin,ymax]
%       Z           - [zmin,zmax]
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
%       P           - simulated wave field
%       x0          - secondary sources
%
%   WAVE_FIELD_MONO_SDM_25D(X,Y,Z,xs,src,f,L,conf) simulates a wave
%   field of the given source type (src) using a SDM 2.5 dimensional driving
%   function in the space/time-frequency domain. This means by calculating 
%   the integral for P with a summation.
%   To plot the result use plot_wavefield(x,y,z,P).
%
%   References:
%       
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: plot_wavefield, wave_field_mono_sdm_25d_kx

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
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
nargmin = 7;
nargmax = 8;
narginchk(nargmin,nargmax);
isargvector(X,Y,Z);
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
% Get the position of the loudspeakers
x0 = secondary_source_positions(L,conf);
% Generate tapering window
win = tapering_window(x0,conf);
% Driving function D(x0,omega)
D = driving_function_mono_sdm_25d(x0,xs,src,f,conf) .* win;
% disable plotting to handle the tapering window
conf.useplot = 0;
% Wave field
[x,y,z,P] = wave_field_mono(X,Y,Z,x0,'ps',D,f,conf);


% ===== Plotting =========================================================
if(useplot)
    plot_wavefield(x,y,z,P,x0,win,conf);
end
