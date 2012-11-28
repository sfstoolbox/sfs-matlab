function label = print_label(dim,unit,conf)
%PRINT_LABEL retruns a suitable axis label
%
%   Usage: label = print_label(dim,[unit,[conf]])
%
%   Input parameters:
%       dim     - name of label dimension
%       unit    - name of label unit, default: no unit
%       conf    - optional configuration struct (see SFS_config)
%
%   Ouput parameters:
%       label   - label to put on an axis of a plot
%
%   PRINT_LABEL(DIM,UNIT) generates a label with the given axis dimension and
%   unit. The formatting is depending on your plotting style, adding $$ for
%   LaTeX labels.
%
%   see also:

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

% FIXME: is this function needed anymore?


%% ===== Checking of input parameter =====================================
nargmin = 1;
nargmax = 3;
narginchk(nargmin,nargmax);
if nargin==nargmax-1
    if isstruct(unit)
        conf = unit;
        unit = '';
    else
        conf = SFS_config;
    end
elseif nargin==nargmax-2
    unit = '';
    conf = SFS_config;
end
isargchar(dim,unit);
isargstruct(conf);


%% ===== Configuration ===================================================
p.mode = conf.plot.mode;


%% ===== Main ============================================================
if strcmp(p.mode,'paper') | strcmp(p.mode,'talk')
    if length(unit)>0
        label = sprintf('$%s$~(%s)',dim,unit);
    else
        label = sprintf('$%s$',dim);
    end
else
    if length(unit)>0
        label = sprintf('%s (%s)',dim,unit);
    else
        label = sprintf('%s',dim);
    end
end
