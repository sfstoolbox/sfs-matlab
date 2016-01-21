function [diam] = secondary_source_maximum_distance(x0)
%SECONDARY_SOURCE_MAXIMUM_DISTANCE calculates the maximum distance 
% between the secondary sources
% (i.e. the diameter of the SSD in Euklidian norm)
%
%   Usage: diam = secondary_source_maximum_distance(x0)
%
%   Input parameters:
%       x0          - secondary sources / m
%
%   Output parameters:
%       diam        - diameter of secondary source distribution / m
%
%   SECONDARAY_SOURCE_MAXIMUM_DISTANCE(x0) calculates the maximum
%   Euklidian distance between the given secondary sources. 
%
%   See also: driving_function_imp_wfs

%*****************************************************************************
% Copyright (c) 2010-2016 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2016 Institut fuer Nachrichtentechnik                   *
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
nargmax = 1;
narginchk(nargmin,nargmax),
isargsecondarysource(x0);


%% ===== Calculation ====================================================
% find the source with largest distance from origin
[~,idx] = max(vector_norm(x0(:,1:3),2));
% find the maximum distace to that source
diam = max(vector_norm(x0(:,1:3) - ...
    repmat(x0(idx,1:3),[size(x0,1),1]),2));
end
