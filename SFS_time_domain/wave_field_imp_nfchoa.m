function varargout = wave_field_imp_nfchoa(X,Y,Z,xs,src,t,conf)
%WAVE_FIELD_IMP_NFCHOA returns the wave field in time domain for NFC-HOA
%
%   Usage: [p,x,y,z,x0] = wave_field_imp_nfchoa(X,Y,Z,xs,src,t,[conf])
%
%   Input options:
%       X           - x-axis / m; single value or [xmin,xmax]
%       Y           - y-axis / m; single value or [ymin,ymax] 
%       Z           - z-axis / m; single value or [zmin,zmax]
%       xs          - position of point source / m
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs, ys are the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%       t           - time point t / samples
%       conf        - optional configuration struct (see SFS_config)
%
%   Output options:
%       p           - simulated wave field
%       x           - corresponding x axis / m
%       y           - corresponding y axis / m
%       z           - corresponding z axis / m
%       x0          - positions and directions of the secondary sources / m
%
%   WAVE_FIELD_IMP_NFCHOA(X,Y,Z,xs,src,t,conf) simulates a wave field of the
%   given source type (src) using a NFC-HOA driving function at the time point t.
%
%   To plot the result use:
%   conf.plot.usedb = 1;
%   plot_wavefield(p,x,y,z,x0,conf);
%
%   see also: driving_function_imp_nfchoa, wave_field_mono_nfchoa

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
isargchar(src);
isargscalar(t);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Computation =====================================================
% Get secondary sources
x0 = secondary_source_positions(conf);
% Calculate driving function
d = driving_function_imp_nfchoa(x0,xs,src,conf);
% Calculate wave field
[varargout{1:min(nargout,4)}] = wave_field_imp(X,Y,Z,x0,'ps',d,t,conf);
% Return secondary sources if desired
if nargout==5 varargout{5}=x0; end
