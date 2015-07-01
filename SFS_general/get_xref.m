function xref = get_xref(x0)
%GET_XREF returns positions of xref positions for 2.5D synthesis
%
%   Usage: xref = get_xref(x0)
%
%   Input parameters:
%       x0      - positions and directions of secondary sources / m [nx6]
%
%   Output parameters:
%       xref    - position of xref points / m [nx3]
%
%   GET_XREF(x0) returns the positions of xref that are used to calculate the 
%   amplitude correction factor g0 for 2.5D synthesis. One position for every
%   secondary source is returned. xref is calculated by shifting the positions
%   of the secondary sources towards their directions.
%
%   See also: driving_function_mono_wfs_pw

%*****************************************************************************
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
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


%% ===== Checking input parameters =======================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
isargmatrix(x0);


%% ===== Main ============================================================
% Create xref as a shifted copy of the secondary source geometry
xref = x0(:,1:3) + x0(:,4:6);
