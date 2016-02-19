function [diam,center] = secondary_source_diameter(conf)
%SECONDARY_SOURCE_DIAMETER calculates the maximum distance
% between the secondary sources (the diameter) and the center of the
% smallest ball that contains the array.
%
%   Usage: [diam,center] = secondary_source_diameter(conf)
%
%   Input parameters:
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       diam        - diameter of secondary source distribution / m
%       center      - center of the ball containing SSD / m [1x3]
%
%   SECONDARAY_SOURCE_DIAMETER(conf) calculates the maximum
%   Euklidian distance between the given secondary sources. Additionaly,
%   the center of the encompassing ball is returned. If one of the predefined
%   secondary source distributions 'linear', 'circular', or 'spherical' is used,
%   the returned diameter is equal to conf.secondary_sources.size.
%   If 'box' or 'rounded-box' is used, the diameter of the bounding box
%   is returned.
%
%   See also: driving_function_imp_wfs, secondary_source_positions

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
isargstruct(conf);


%% ===== Configuration ==================================================
geometry = conf.secondary_sources.geometry;


%% ===== Calculation ====================================================
if strcmp('line',geometry)   || strcmp('linear',geometry)    || ...
   strcmp('circle',geometry) || strcmp('circular',geometry)  || ...
   strcmp('sphere',geometry) || strcmp('spherical',geometry)
    diam = conf.secondary_sources.size;
    center = conf.secondary_sources.center;
elseif strcmp('box',geometry)
    dx0 = conf.secondary_sources.size/(conf.secondary_sources.number/4-1);
    diam = (conf.secondary_sources.size+dx0)*sqrt(2);
    center = conf.secondary_sources.center;
elseif strcmp('rounded-box',geometry)
    diam = (conf.secondary_sources.size)*sqrt(2);
    center = conf.secondary_sources.center;
else
    x0 = conf.secondary_sources.x0;
    if isempty(x0)
        error(['%s: conf.secondary_sources.x0 must contain the secondary ',...
            'sources when using geometry %s.'],upper(mfilename),geometry);
    end
    % Find source1 :=  source with largest distance from origin
    [~,idx1] = max(vector_norm(x0(:,1:3),2));
    % Find source2 := source with maximum distace to source1
    [diam,idx2] = max(vector_norm(x0(:,1:3) - ...
        repmat(x0(idx1,1:3),[size(x0,1),1]),2));
    % Center is half-way between source1 and source2
    center = x0(idx1,1:3) +  0.5 * (x0(idx2,1:3) - x0(idx1,1:3));
end
