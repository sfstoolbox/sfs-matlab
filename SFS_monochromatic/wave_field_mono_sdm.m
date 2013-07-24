function varargout = wave_field_mono_sdm(X,Y,Z,xs,src,f,conf)
%WAVE_FIELD_MONO_SDM simulates a wave field for WFS
%
%   Usage: [P,x,y,z,x0] = wave_field_mono_sdm(X,Y,Z,xs,src,f,[conf])
%
%   Input parameters:
%       X           - x-axis / m; single value or [xmin,xmax]
%       Y           - y-axis / m; single value or [ymin,ymax]
%       Z           - z-axis / m; single value or [zmin,zmax]
%       xs          - position of virtual source / m
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs is the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       f           - monochromatic frequency / Hz
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       P           - simulated wave field
%       x           - corresponding x axis / m
%       y           - corresponding y axis / m
%       z           - corresponding z axis / m
%       x0          - active secondary sources / m
%
%   WAVE_FIELD_MONO_SDM(X,Y,Z,xs,src,f,conf) simulates a wave field for the
%   given source type (src) using SDM driving functions in the temporal domain.
%   This means by calculating the integral for P with a summation.
%   To plot the result use plot_wavefield(P,x,y,z,x0,win).
%
%   References:
%       Spors2009 - Physical and Perceptual Properties of Focused Sources in
%           Wave Field Synthesis (AES127)
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: plot_wavefield, wave_field_imp_wfs, driving_function_mono_wfs

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
nargmin = 6;
nargmax = 7;
narginchk(nargmin,nargmax);
isargvector(X,Y,Z);
isargxs(xs);
isargpositivescalar(f);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
if strcmp('2D',conf.dimension)
    greens_function = 'ls';
else
    greens_function = 'ps';
end


%% ===== Computation ====================================================
% Get the position of the loudspeakers and its activity
x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,xs,src);
x0 = secondary_source_tapering(x0,conf);
% Driving function
D = driving_function_mono_sdm(x0,xs,src,f,conf);
% Wave field
[varargout{1:min(nargout,4)}] = ...
    wave_field_mono(X,Y,Z,x0,greens_function,D,f,conf);
% Return secondary sources if desired
if nargout==5, varargout{5}=x0; end
