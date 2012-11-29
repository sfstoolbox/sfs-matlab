function addirspath(varargin)
%ADDIRSPATH adds directories containing irs files to the path
%
%   Usage: addirspath(varargin)
%
%   Input parameters:
%       varargin - path or pathe containing irs data sets.
%                  Default: '~/svn/ir_databases' and '~/svn/measurements'
%
%   ADDIRSPATH(varargin) adds the given directorysand its subdirectories to
%   the path. If no directory is given, '~/svn/ir_databases' and
%   '~/svn/measurements' are added.

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


% ===== Checking of input parameters ====================================
nargmin = 0;
nargmax = inf;
narginchk(nargmin,nargmax);
if nargin==nargmin
    dirs{1} = '~/svn/ir_databases';
    dirs{2} = '~/svn/measurements';
else
    isargdir(char(varargin));
    dirs = varargin;
end


%% ===== Adding pathes ==================================================
for ii = 1:length(dirs)
    path = genpath(dirs{ii});
    addpath(path,'-end');
end
