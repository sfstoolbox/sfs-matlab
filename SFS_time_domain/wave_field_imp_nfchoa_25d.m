function varargout = wave_field_imp_nfchoa_25d(X,Y,Z,xs,src,t,L,conf)
%WAVE_FIELD_IMP_NFCHOA_25D returns the wave field in time domain of an impulse
%
%   Usage: [p,x,y,z,x0] = wave_field_imp_nfchoa_25d(X,Y,Z,xs,src,t,L,[conf])
%
%   Input options:
%       X           - [xmin,xmax]
%       Y           - [ymin,ymax]
%       Z           - [zmin,zmax]
%       xs          - position of point source (m)
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs, ys are the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%       t           - time point t (samples)
%       L           - array length (m)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output options:
%       p           - simulated wave field
%       x           - corresponding x axis
%       y           - corresponding y axis
%       z           - corresponding z axis
%       x0          - positions and directions of the secondary sources
%
%   WAVE_FIELD_IMP_NFCHOA_25D(X,Y,Z,xs,src,t,L,conf) simulates a wave field of the
%   given source type (src) using a NFC-HOA 2.5 dimensional driving
%   function at the time point t.
%
%   To plot the result use:
%   x0 = secondary_source_positions(L,conf);
%   conf.plot.usedb = 1;
%   plot_wavefield(p,x,y,z,x0,conf);

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
isargxs(xs);
isargpositivescalar(L);
isargchar(src);
isargscalar(t);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Computation =====================================================
% Get secondary sources
x0 = secondary_source_positions(L,conf);
% Calculate driving function
d = driving_function_imp_nfchoa_25d(x0,xs,src,L,conf);
% Calculate wave field
[varargout{1:min(nargout,4)}] = wave_field_imp_3d(X,Y,Z,x0,d,t,conf);
% Return secondary sources if desired
if nargout==5 varargout{5}=x0; end
