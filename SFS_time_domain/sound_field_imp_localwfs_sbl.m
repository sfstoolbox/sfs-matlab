function varargout = sound_field_imp_localwfs_sbl(X,Y,Z,xs,src,t,conf)
%SOUND_FIELD_IMP_LOCALWFS_SBL returns the sound field in time domain of an 
%impulse
%
%   Usage: [p,x,y,z,x0] = sound_field_imp_localwfs_sbl(X,Y,Z,xs,src,t,conf)
%
%   Input options:
%       X           - x-axis / m; single value or [xmin,xmax] or nD-array
%       Y           - y-axis / m; single value or [ymin,ymax] or nD-array
%       Z           - z-axis / m; single value or [zmin,zmax] or nD-array
%       xs          - position of point source / m
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs, ys are the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%       t           - time point t of the sound field / seconds
%       conf        - configuration struct (see SFS_config)
%
%   Output options:
%       p           - simulated sound field
%       x           - corresponding x values / m
%       y           - corresponding y values / m
%       z           - corresponding z values / m
%       x0          - secondary sources / m
%
%   SOUND_FIELD_IMP_LOCALWFS(X,Y,Z,xs,src,t,conf) simulates a sound field of the
%   given source type (src) synthesized with local wave field synthesis at the
%   time t.
%
%   To plot the result use:
%   plot_sound_field(p,X,Y,Z,x0,conf);
%   or simple call the function without output argument:
%   sound_field_imp_localwfs_sbl(X,Y,Z,xs,src,t,conf)
%   For plotting you may also consider to display the result in dB, by setting
%   the following configuration option before:
%   conf.plot.usedB = true;
%
%   See also: driving_function_imp_localwfs_sbl

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


%% ===== Checking of input  parameters ==================================
nargmin = 7;
nargmax = 7;
narginchk(nargmin,nargmax);
isargnumeric(X,Y,Z);
isargxs(xs);
isargchar(src);
isargscalar(t);
isargstruct(conf);
fs = conf.fs;


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
[d,delay_offset] = driving_function_imp_localwfs_sbl(x0,xs,src,conf);
% Ensure virtual source/secondary source activity starts at t = 0
t = t + delay_offset;
% Calculate sound field
[varargout{1:min(nargout,4)}] = ...
    sound_field_imp(X,Y,Z,x0,greens_function,d,t,conf);
% Return secondary sources if desired
if nargout==5, varargout{5}=x0; end
