function [number,L] = secondary_source_number(L,dx0,conf)
%SECONDARY_SOURCE_NUMBER calculate the number of secondary sources for the
%   given loudspeaker distance
%
%   Usage: number = secondary_source_number(dx0,[conf])
%
%   Input parameters:
%       L       - length/diameter of secondary source array / m
%       dx0     - distance between secondary sources / m
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       number  - number of needed loudspeaker
%       L       - corrected array size, corresponding to the given dx0
%
%   SECONDARAY_SOURCE_NUMBER(L,dx0,conf) calculates the number of secondary
%   sources for the given distance dx0 and array size, considering the array
%   geometry conf.secondary_sources.geometry. It returns also a corrected array
%   size L, because for some geometries like a circle it could be that the given
%   distance dx0 does not fit completely to the given perimeter. 
%
%   see also: secondary_source_positions, secondary_source_selection

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 3;
narginchk(nargmin,nargmax),
isargpositivescalar(dx0,L);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
% Array type
geometry = conf.secondary_sources.geometry;
% Predefined secondary sources
x0 = conf.secondary_sources.x0;


%% ===== Calculation ====================================================
% Check if we have given secondary sources
if ~isempty(x0)
    isargsecondarysource(x0);
    number = size(x0,1);
elseif strcmp('linear',geometry)
    % Number of secondary sources
    number = fix(L/dx0)+1;
    % Corresponding size of loudspeaker array
    L = (number-1)*dx0;
elseif strcmp('circle',geometry) || strcmp('circular',geometry) || ...
       strcmp('spherical',geometry) || strcmp('sphere',geometry)
    % L is the diameter!
    % Perimeter of the circle
    P = pi*L;
    % Number of loudspeakers
    number = round(P/dx0);
    % Corresponding size of loudspeaker array
    L = (number*dx0)/pi;
elseif strcmp('box',geometry)
    % Number of loudspeakers
    number = 4*(fix(L/dx0)+1);
    % Corresponding size of loudspeaker array
    L = (number/4-1)*dx0;
else
    error('%s: %s is a unknown array geometry.',upper(mfilename),geometry);
end
