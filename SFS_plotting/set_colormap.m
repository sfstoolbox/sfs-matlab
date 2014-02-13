function set_colormap(map)
%SET_COLORMAP set the color of the color map
%
%   Usage: set_colormap(map)
%
%   Input options:
%       map         - can be the real map (see help colormap) or one of the
%                     following:
%                       'default' - blue,white,red
%                       'gray'    - white,black 
%
%   SET_COLORMAP(map) sets the color of the color map to map. map can be a
%   three column matrix containing the R,G,B values between 0 and 1 or it can
%   be a string to choose from a defined set of available maps.
%
%   see also: set_colorbar, plot_wavefield

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


%% ===== Checking of input parameter =====================================
nargmin = 0;
nargmax = 1;
narginchk(nargmin,nargmax);
if nargin<nargmax
    map = 'default';
end


%% ===== Plotting ========================================================
% Change color map (default: gray)
if ~ischar(map)
    colormap(map);
elseif strcmp('default',map) || strcmp('moreland',map)
    % Default SFS Toolbox color: blue,white,red
    % see: http://www.sandia.gov/~kmorel/documents/ColorMaps/
    colormap(moreland(256));
elseif strcmp('gray',map) || strcmp('grey',map) 
    colormap(flipud(colormap('gray')));
else
    colormap(map);
end
