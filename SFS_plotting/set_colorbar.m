function set_colorbar(conf)
%SET_COLORBAR draws a color bar to the plot
%
%   Usage: set_colorbar([conf])
%
%   Input options:
%       conf        - optional configuration struct (see SFS_config)
%
%   SET_COLORBAR() drwas a color bar on the figure and sets the map to the color
%   specified in conf.plot.colormap.
%
%   see also: plot_wavefield, set_colormap

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
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
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
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
p.caxis = conf.plot.caxis;
p.usedb = conf.plot.usedb;
p.colormap = conf.plot.colormap;


%% ===== Plotting ========================================================
% Change color map (default: gray)
set_colormap(p.colormap);

% Set the limits of the colormap and add a colorbar
if length(p.caxis)==2
    caxis(p.caxis);
end
if p.usedb
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
end


