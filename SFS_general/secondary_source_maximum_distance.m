function [diam,center] = secondary_source_maximum_distance(x0)
%SECONDARY_SOURCE_MAXIMUM_DISTANCE calculates the maximum distance 
% between the secondary sources (the diameter) and the center of the 
% smallest ball that contains the array.
% 
%
%   Usage: [diam,center] = secondary_source_maximum_distance(x0)
%
%   Input parameters:
%       x0          - secondary sources / m [nx7]
%
%   Output parameters:
%       diam        - diameter of secondary source distribution / m
%       center      - center of the ball containing SSD / m [1x3]   
%   
%   SECONDARAY_SOURCE_MAXIMUM_DISTANCE(x0) calculates the maximum
%   Euklidian distance between the given secondary sources. Additionaly,
%   the center of the encompassing is returned.
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
narginchk(nargmin,nargmax);

%% ===== Calculation ====================================================
% find source1 :=  source with largest distance from origin
[~,idx1] = max(vector_norm(x0(:,1:3),2));
% find source2 := source with maximum distace to source1
[diam,idx2] = max(vector_norm(x0(:,1:3) - ...
    repmat(x0(idx1,1:3),[size(x0,1),1]),2));
% center is half-way between source1 and source2
center = x0(idx1,1:3) +  0.5 * (x0(idx2,1:3) - x0(idx1,1:3));
end
