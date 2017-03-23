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
%   See also: set_colorbar, plot_sound_field

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
