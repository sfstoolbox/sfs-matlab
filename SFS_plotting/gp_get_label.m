function label = gp_get_label(dim,unit,terminal)
%GP_GET_LABEL returns a suitable axis label for gnuplot
%
%   Usage: label = gp_get_label(dim,[unit,[terminal]])
%
%   Input parameters:
%       dim         - name of label dimension
%       unit        - name of label unit (default: '')
%       terminal    - name of the terminal (default: 'wxt')
%                     available terminals are: 'wxt','png','eps','epslatex'
%
%   Ouput parameters:
%       label   - label to put on an axis of a plot
%
%   GP_GET_LABEL(dim,unit,terminal) generates a label with the given axis
%   dimension and unit. The formatting is depending on the given terminal, 
%   adding $$ for LaTeX labels.
%
%   see also: gp_print_screen, gp_print_png, gp_print_eps, gp_print_epslatex

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
nargmin = 1;
nargmax = 3;
narginchk(nargmin,nargmax);
if nargin==nargmax-1
    terminal = 'wxt';
elseif nargin==nargmax-2
    unit = '';
    terminal = 'wxt';
end
isargchar(dim,unit,terminal);


%% ===== Main ============================================================
if strcmp('epsterminal',terminal)
    if length(unit)>0
        label = sprintf('$%s /$\;%s',dim,unit);
    else
        label = sprintf('$%s$',dim);
    end
else
    if length(unit)>0
        label = sprintf('%s / %s',dim,unit);
    else
        label = sprintf('%s',dim);
    end
end
