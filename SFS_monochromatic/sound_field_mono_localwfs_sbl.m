function varargout = sound_field_mono_localwfs_sbl(X,Y,Z,xs,src,f,conf)
%SOUND_FIELD_MONO_LOCALWFS_SBL sound field of local WFS using spatial bandwith
%limitation
%
%   Usage: [P,x,y,z,x0] = sound_field_mono_localwfs_sbl(X,Y,Z,xs,src,f,conf)
%
%   Input options:
%       X           - x-axis / m; single value or [xmin,xmax] or nD-array
%       Y           - y-axis / m; single value or [ymin,ymax] or nD-array
%       Z           - z-axis / m; single value or [zmin,zmax] or nD-array
%       xs          - position of point source or direction of plane
%                     wave / m [1x3]
%       src         - source type of the virtual source
%                         'pw' - plane wave
%                         'ps' - point source
%       f           - frequency of the monochromatic source / Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output options:
%       P           - simulated sound field
%       x           - corresponding x values / m
%       y           - corresponding y values / m
%       z           - corresponding z values / m
%       x0          - secondary sources / m
%
%   SOUND_FIELD_MONO_LOCALWFS_SBL(X,Y,Z,xs,src,f,conf) simulates a monochromatic
%   sound field of the given source type (src) synthesized with local wave field
%   synthesis using spatial bandwidth limitation for the frequency f.
%
%   To plot the result use:
%   plot_sound_field(P,X,Y,Z,x0,conf);
%   or simple call the function without output argument:
%   sound_field_mono_localwfs_sbl(X,Y,Z,xs,src,f,conf)
%   For plotting you may also consider to display the result in dB, by setting
%   the following configuration option before:
%   conf.plot.usedB = true;
%
%   See also: plot_sound_field, driving_function_mono_localwfs_sbl

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2019 SFS Toolbox Developers                             *
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
% https://sfs.readthedocs.io                            sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 7;
nargmax = 7;
narginchk(nargmin,nargmax);
isargnumeric(X,Y,Z);
isargxs(xs);
isargchar(src);
isargscalar(f);
isargstruct(conf);


%% ===== Configuration ==================================================
if strcmp('2D',conf.dimension)
    greens_function = 'ls';
else
    greens_function = 'ps';
end


%% ===== Computation =====================================================
% Get secondary sources
x0 = secondary_source_positions(conf);
% Get driving signals
D = driving_function_mono_localwfs_sbl(x0,xs,src,f,conf);
% Calculate sound field
[varargout{1:min(nargout,4)}] = ...
    sound_field_mono(X,Y,Z,x0,greens_function,D,f,conf);
% Return secondary sources if desired
if nargout==5, varargout{5}=x0; end
