function x0 = secondary_source_tapering(x0,varargin)
%SECONDARY_SOURCE_TAPERING applies a tapering window to the secondary sources
%
%   Usage: x0 = secondary_source_tapering(x0,conf)
%
%   Input options:
%       x0          - secondary sources / m
%       conf        - optional configuration struct (see SFS_config)
%
%   Output options:
%       x0          - secondary sources / m, containing the applied tapering
%                     window in its weights in x0(:,7)
%
%   SECONDARY_SOURCE_TAPERING(x0,conf) applies a tapering window (depending on
%   the conf.usetapwin and conf.tapwinlen settings) to the secondary sources.
%   It is applied to the weights stored in x0(:,7).
%
% see also: secondary_source_positions, tapering_window

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
narginchk(nargmin,nargmax);


%% ===== Calculation ====================================================
% Apply tapering window to secondary sources
x0(:,7) = x0(:,7) .* tapering_window(x0,varargin{:});
