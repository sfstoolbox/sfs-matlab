function bool = is_dim_custom(varargin)
%IS_DIM_CUSTOM returns true for a custom grid, otherwise false
%
%   Usage: bool = is_grid_custom(x1,x2,...)
%
%   Input parameters:
%       x1,x2,... - axis / m; single value or [xmin,xmax] or nD-array
%
%   Output parameters:
%       bool      - array of logical indicating whether each input is an 
%                   nD-array
%
%   IS_DIM_CUSTOM(x1,x2,..) checks if we have a custom grid by checking if any 
%   of the given x,y,z values is a nD-array.
%
%   See also: xyz_grid, xyz_axes_selection, plot_sound_field

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
isargnumeric(varargin{:});


%% ===== Computation =====================================================
bool = cellfun(@(x) numel(x)>2, varargin);
