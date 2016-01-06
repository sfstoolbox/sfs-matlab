function map = chromajs(n)
%CHROMAJS returns a divergent lightyellow-red color map
%
%   Usage: map = chromajs([n])
%
%   Input parameters:
%       n   - optional length of colormap (default uses the figure default
%             length)
%
%   Output parameters:
%       map - colormap [n 3]
%
%   CHROMAJS(N) returns an N-by-3 matrix containing a divergent colormap.
%   For details on the colormap have a look at:
%   http://gka.github.io/palettes/#colors=lightyellow,orangered,deeppink,darkred|steps=7|bez=1|coL=1
%
%   For example, to reset the colormap of the current figure:
%             colormap(chromajs)
%
% See also: generate_colormap, moreland, plot_sound_field

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


%% ===== Checking input parameters =======================================
nargmin = 0;
nargmax = 1;
narginchk(nargmin,nargmax);
if nargin < nargmax
    n = size(get(gcf,'colormap'),1);
end


%% ===== Computation =====================================================
table = [ 255, 255, 224
255, 223, 184
255, 188, 148
255, 151, 119
255, 105,  98
238,  66,  86
210,  31,  71
176,   6,  44
139,   0,   0];
map = generate_colormap(table,n);
