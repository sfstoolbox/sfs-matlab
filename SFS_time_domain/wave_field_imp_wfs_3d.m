function [x,y,z,p,x0,win] = wave_field_imp_wfs_3d(X,Y,Z,xs,src,t,L,conf)
%WAVE_FIELD_IMP_WFS_3D returns the wave field in time domain of an impulse
%
%   Usage: [x,y,p,dds] = wave_field_imp_wfs_3d(X,Y,xs,src,t,L,[conf])
%
%   Input options:
%       X           - length of the X axis (m); single value or [xmin,xmax]
%       Y           - length of the Y axis (m); single value or [ymin,ymax]
%       Z           - length of the Z axis (m); single value or [zmin,zmax]
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
%       x,y,z       - x-,y- and z-axis of the wave field
%       p           - wave field (length(y) x length(x))
%       x0          - active secondary sources
%       win         - tapering window
%
%   WAVE_FIELD_IMP_WFS_3D(X,Y,xs,src,t,L,conf) simulates a wave field of the
%   given source type (src) using a WFS 3 dimensional driving function with
%   a delay line.
%   To plot the result use:
%   conf.plot.usedb = 1;
%   plot_wavefield(x,y,z,p,x0,win,conf);

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

% AUTHOR: Hagen Wierstorf, Sascha Spors
% $LastChangedDate: 2012-06-14 11:03:36 +0200 (Thu, 14 Jun 2012) $
% $LastChangedRevision: 761 $
% $LastChangedBy: wierstorf.hagen $


%% ===== Checking of input  parameters ==================================
nargmin = 7;
nargmax = 8;
error(nargchk(nargmin,nargmax,nargin));
isargvector(X,Y,Z);
xs = position_vector(xs);
isargpositivescalar(L);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end

%% ===== Configuration ==================================================
xref = position_vector(conf.xref);
useplot = conf.useplot;


%% ===== Computation =====================================================
% Get secondary sources
x0 = secondary_source_positions(L,conf);
x0 = secondary_source_selection(x0,xs,src,xref);
% Generate tapering window
win = tapering_window(x0,conf);

% Get driving signals
d = driving_function_imp_wfs_3d(x0,xs,src,conf);
% Apply tapering window
d = bsxfun(@times,d,win');

% disable plotting in order to integrate the tapering window
conf.useplot = 0;
% Calculate wave field
[x,y,z,p] = wave_field_imp_3d(X,Y,Z,x0,d,t,conf);


%% ===== Plotting ========================================================
if useplot
    conf.plot.usedb = 1;
    plot_wavefield(x,y,z,p,x0,win,conf);
end
