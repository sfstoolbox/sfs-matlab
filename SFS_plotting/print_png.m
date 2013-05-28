function print_png(outfile,conf);
%FIGSIZE changes the size of a figure
%
%   Usage: figsize(x,y,unit)
%
%   Input options:
%       x,y         - x,y size of the figure
%       unit        - unit in which the size is given, can be one of the
%                     following: 'cm', 'px'
%
%   FIGSIZE(x,y,unit) sets the size of the last figure to x,y in the given unit.
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
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
p.resolution = conf.plot.resolution;
p.size = conf.plot.size;
p.size_unit = conf.plot.size_unit;


%% ===== Main ============================================================
dpi = sprintf('-r%i',p.resolution);
if isoctave
    if ~strcmp('px',p.size_unit)
        error('%s: unit has to be in px under Octave for a png plot', ...
            upper(mfilename));
    end
    res = sprintf('-S%i,%i',p.size(1),p.size(2));
    print(outfile,'-dpng',dpi,res);
else
    print(outfile,'-dpng',dpi);
end
close;

