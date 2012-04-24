function print_png(filename)
%PRINT_PNG creates a png file from the last plot command
%
%   Usage: print_png(filename)
%
%   Input options:
%       filename    - name of the output png file
%
%   PRINT_PNG(filename) creates a png version from the last plot command in the
%   in the given file. The size of the png will be 500x375 pixel. The line width
%   etc. of the axis is arranged in order to create a nice looking plot.
%
%   see also: print_eps

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

% Is this function needed anymore?

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
if ~ischar(filename)
    error('%s: filename has to be a string!',upper(mfilename));
end


%% ===== Plotting ========================================================
% Save old values in order to restore them
linewidth = get(gca,'LineWidth');
ticklength = get(gca,'TickLength');
position = get(gca,'Position');
set(gca,'LineWidth',2);
set(gca,'TickLength',[0.005 0.005]);
% FIXME: check if the following also doesn't work under Matlab!
% Update Position by hand to aviod clipping of labels
%sc = 1.2;
%set(gca,'Position',...
%    [position(1)*sc,position(2)*sc,position(3)/sc,position(4)/sc]);
print(filename,'-dpng','-r75');
set(gca,'LineWidth',linewidth);
set(gca,'TickLength',ticklength);
set(gca,'Position',position);
