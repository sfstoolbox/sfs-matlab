function [p,x,y,z,x0,win] = wave_field_imp_wfs(X,Y,Z,xs,src,t,L,conf)
%WAVE_FIELD_IMP_WFS returns the wave field in time domain of an impulse
%
%   Usage: [p,x,y,z,x0,win] = wave_field_imp_wfs(X,Y,Z,xs,src,t,L,[conf])
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
%                         'fs' - focused source
%       t           - time point t of the wave field (samples)
%       L           - array length (m)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output options:
%       p           - simulated wave field
%       x           - corresponding x axis
%       y           - corresponding y axis
%       z           - corresponding z axis
%       x0          - secondary sources
%       win         - tapering window
%
%   WAVE_FIELD_IMP_WFS(X,Y,Z,xs,src,t,L,conf) simulates a wave field of the
%   given source type (src) using a WFS driving function with a delay line at
%   the time t.
%
%   To plot the result use:
%   conf.plot.usedb = 1;
%   plot_wavefield(p,x,y,z,x0,win,conf);

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


%% ===== Configuration ==================================================
xref = conf.xref;
useplot = conf.useplot;


%% ===== Computation =====================================================
% Get secondary sources
x0 = secondary_source_positions(L,conf);
x0 = secondary_source_selection(x0,xs,src);
% Generate tapering window
win = tapering_window(x0,conf);

% Get driving signals
d = driving_function_imp_wfs(x0,xs,src,conf);
% Apply tapering window
d = bsxfun(@times,d,win');

% disable plotting in order to integrate the tapering window
conf.useplot = 0;
% Calculate wave field
[p,x,y,z] = wave_field_imp(X,Y,Z,x0,'ps',d,t,conf);


%% ===== Plotting ========================================================
if nargout==0 || useplot
    conf.plot.usedb = 1;
    plot_wavefield(p,x,y,z,x0,win,conf);
end
