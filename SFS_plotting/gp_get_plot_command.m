function cmd = gp_get_plot_command(p)
%GP_GET_PLOT_COMMAND returns the gnuplot command for plotting the sound field
%and loudspeaker
%
%   Usage: cmd = gp_get_plot_command(p)
%
%   Input parameters:
%       p   - struct holding all the conf.plot data and more
%
%   Output parameters:
%       cmd - plot command to insert in gnuplot file
%
%   GP_GET_PLOT_COMMAND(p) returns the plot command for gnuplot to plot the
%   sound field and loudspeaker as symbols or points, or completly without
%   loudspeakers.
%
%   See also: gp_print_*, plot_wavefield

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


if p.loudspeakers && p.realloudspeakers
    % Plotting real loudspeaker symbols
    cmd = sprintf([ ...
        'call ''gp_draw_loudspeakers.gnu'' ''%s'' ''%f''\n', ...
        'plot ''%s'' binary matrix with image'], ...
        p.lsfile,p.lssize,p.datafile);
elseif p.loudspeakers
    % Plotting only points at the loudspeaker positions
    cmd = sprintf([ ...
        'plot ''%s'' binary matrix with image,\\\n', ...
        '     ''%s'' u 1:2 w points ls 1'], ...
        p.datafile,p.lsfile);
else
    % Plotting no loudspeakers at all
    cmd = sprintf( ...
        'plot ''%s'' binary matrix with image', ...
        p.datafile);
end
