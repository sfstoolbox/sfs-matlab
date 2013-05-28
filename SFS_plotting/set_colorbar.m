function set_colorbar(conf)
%SET_COLORBAR draws a color bar to the plot
%
%   Usage: draw_loudspeakers([conf])
%
%   Input options:
%       x0          - positions and directions of the loudspeakers / m
%       win         - tapering window, which is the activity of the loudspeaker
%                     (default: 1)
%       conf        - optional configuration struct (see SFS_config)
%
%   DRAW_LOUDSPEAKERS(x0,win) draws loudspeaker symbols or crosses at the given
%   secondary source positions. This can be controlled by the
%   conf.plot.realloudspeakers setting. The loudspeaker symbols are pointing in
%   their given direction. In addition to the secondary source positions, the
%   activity of the single secondary sources can be given by the win vector. For
%   every secondary source it can contain a value between 0 and 1. 1 is fully
%   active. If only one value if given, it is used for all secondary sources.
%
%   see also: plot_wavefield

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


%% ===== Checking of input parameter =====================================
nargmin = 0;
nargmax = 1;
narginchk(nargmin,nargmax);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
p.caxis = plot.caxis;
p.usedb = plot.usedb;


%% ===== Plotting ========================================================
% Set the limits of the colormap and add a colorbar
if length(p.caxis)==2
    caxis(p.caxis);
end
h = colorbar;
ylabel(h,'Amplitude (dB)');
% Get the font size and name of the figure and adjust the colorbar
fsize = get(gca,'FontSize');
fname = get(gca,'FontName');
set(h,'FontSize',fsize);
set(h,'FontName',fname);
temp = get(h,'Ylabel');
set(temp,'FontSize',fsize);
set(temp,'FontName',fname);
