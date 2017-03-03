function set_colorbar(conf)
%SET_COLORBAR draws a color bar to the plot
%
%   Usage: set_colorbar(conf)
%
%   Input options:
%       conf        - configuration struct (see SFS_config)
%
%   SET_COLORBAR(conf) draws a color bar on the figure and sets the map to the
%   color specified in conf.plot.colormap.
%
%   See also: plot_sound_field, set_colormap

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input parameter =====================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
isargstruct(conf);


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


