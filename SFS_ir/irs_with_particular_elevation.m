function irs = irs_with_particular_elevation(irs,delta)
%IRS_WITH_PARTICULAR_ELEVATION(irs,delta) returns an IRS set which only contains
%data for one elevation
%
%   Usage: irs = irs_with_particular_elevation(irs,[delta])
%
%   Input parameters:
%       irs     - IR data set
%       delta   - elevation angle for the desired IR / rad
%                 default: 0
%
%   Output parameters:
%       irs      - IRS for the given elevation
%
%   IRS_WITH_PARTICULAR_ELEVATION(irs,delta) returns a IRS set for the given
%   angle delta, or by default the horizontal plane. Input should be an IRS-set
%   with diffrent values for the elevation
%
%   see also: slice_irs, new_irs

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


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax)
check_irs(irs);
if nargin==nargmax
    isargscalar(delta);
else
    delta = 0;
end


%% ===== Computation ====================================================

% Finding the entries belonging to delta and slice the irs
idx = (( round(irs.apparent_elevation*10)==round(10*delta) ));
irs = slice_irs(irs,idx);
