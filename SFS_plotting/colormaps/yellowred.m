function map = yellowred(n)
%CHROMAJS returns a divergent lightyellow-red color map
%
%   Usage: m = yellowred([n])
%
%   Input parameters:
%       n - optional length of colormap (default uses the figure default length)
%
%   Output parameters:
%       m - colormap [n 3]
%
%   YELLOWRED(N) returns an N-by-3 matrix containing a divergent colormap.
%   Without a given N the same length as the current figure's colormap is used.
%   For details on the colormap have a look at: http://bit.ly/2vC3Ogr
%
%   To change the colormap current figure run: colormap(yellowred)
%
% See also: generate_colormap, moreland, plot_sound_field

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
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
