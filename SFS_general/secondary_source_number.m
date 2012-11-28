function [nls,L] = secondary_source_number(L,conf)
%SECONDARY_SOURCE_NUMBER calculate the number of secondary sources for the
%   given array size
%
%   Usage: [nls,L] = secondary_source_number(L,[conf])
%
%   Input parameters:
%       L       - length of the loudspeaker array (m)
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       nls     - number of needed loudspeaker
%       L       - real length of the loudspeaker array (corresponding to
%                 conf.dx0)
%
%   SECONDARAY_SOURCE_NUMBER(L,conf) calculates the number of needed loudspeaker for
%   the given array length L, using the config loudspeaker distance conf.dx0.
%   Also the real length L of such a loudspeaker array will returned. This is
%   neccessary because the user given length L is probably not compatible to the
%   given value conf.dx0.
%
%   see also: secondary_source_positions, secondary_source_selection,
%       secondary_source_selection

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


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax),
isargpositivescalar(L);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
% Array type
array = conf.array;
% Distance between secondary sources
dx0 = conf.dx0;
% Predefined secondary sources
x0 = conf.x0;


%% ===== Calculation ====================================================
% Check if we have given secondary sources
if ~isempty(x0)
    isargsecondarysource(x0);
    nls = size(x0,1);
elseif strcmp('linear',array)
    % Number of loudspeaker
    nls = fix(L/dx0)+1;
    % Corresponding size of loudspeaker array
    L = (nls-1)*dx0;
elseif strcmp('circle',array)
    % L is the diameter!
    % Perimeter of the circle
    P = pi*L;
    % Number of loudspeakers
    nls = round(P/dx0);
    % Corresponding size of loudspeaker array
    L = (nls*dx0)/pi;
elseif strcmp('box',array)
    % Number of loudspeakers
    nls = 4*(fix(L/dx0)+1);
    % Corresponding size of loudspeaker array
    L = (nls/4-1)*dx0;
else
    error('%s: %s is a unknown array type.',upper(mfilename),array);
end
