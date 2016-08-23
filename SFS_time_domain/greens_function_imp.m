function [g,t] = greens_function_imp(x,y,z,xs,src,t,conf)
%GREENS_FUNCTION_IMP returns a Green's function in the time domain
%
%   Usage: [g,t] = greens_function_imp(x,y,z,xs,src,t,conf)
%
%   Input parameters:
%       x       - x points / m
%       y       - y points / m
%       z       - z points / m
%       xs      - position of the source / m
%       src     - source model of the Green's function. Valid models are:
%                   'ps' - point source
%                   'ls' - line source
%                   'pw' - plane wave
%       t       - time / s
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       g       - Green's function evaluated at the points x,y,z
%       t       - Corresponding time values with integrated time
%                 shift / s
%
%   GREENS_FUNCTION_IMP(x,y,z,xs,src,t,conf) calculates the Green's function for
%   the given source model located at xs for the given points x,y,z. Note, that
%   the delta function for the time t is returned as an extra argument. If you
%   want the value of the Green's function only for a specific time you should
%   have a look at sound_field_imp() and apply the folowing command:
%   [p,x,y,z] = sound_field_imp(X,Y,Z,[xs 0 -1 0],src,1,t,conf);
%
%   See also: greens_function_mono, sound_field_imp

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
% disabled checking for performance reasons


%% ===== Configuration ===================================================
c = conf.c;


%% ===== Computation =====================================================
% calculate Green's function for the given source model
if strcmp('ps',src)
    % Source model for a point source: 3D Green's function.
    %
    %                  1
    % g(x-xs,t) = ---------- delta(t - |x-xs|/c)
    %             4pi |x-xs|
    %
    % See http://sfstoolbox.org/#equation-s.ps
    %
    r = sqrt((x-xs(1)).^2+(y-xs(2)).^2+(z-xs(3)).^2);
    g = 1./(4*pi.*r);
    t = (r/c)-t;

elseif strcmp('dps',src)
    % Source model for a dipole point source: derivative of 3D Green's function.
    %
    %                 1   / -1 / iw \      1    \  (x-xs)ns
    % g(x-xs,ns,t) = --- | F  | ---- | + ------  | --------- delta(t - |x-xs|/c)
    %                4pi  \    \ c  /    |x-xs| /  |x-xs|^2
    %
    % See http://sfstoolbox.org/#equation-s.dps
    %
    to_be_implemented(mfilename);

elseif strcmp('ls',src)
    % Source model for a line source: 2D Green's function.
    %                          ___
    %              -1/ c\     | 1       1
    % g(x-xs,t) = F |--  |  - |---  --_-_-_- delta(t - |x-xs|/c)
    %                \iw/    \|8pi  \||x-xs|
    %
    % See http://sfstoolbox.org/en/latest/#equation-s.ls
    % Note, that the filter F^-1 is not implemented!!!!
    %
    r = sqrt((x-xs(1)).^2+(y-xs(2)).^2+(z-xs(3)).^2);
    g = 1./sqrt(r) * sqrt(1/(8*pi));
    t = (r/c)-t;

elseif strcmp('pw',src)
    % Source model for a plane wave:
    %
    % g(x,t) = delta(t - nx/c)
    %
    % See http://sfstoolbox.org/#equation-s.pw
    %
    % direction of plane wave
    nxs = xs / norm(xs);
    %
    g = 1;
    t = (nxs(1).*x+nxs(2).*y+nxs(3).*z)./c-t;

else
    error('%s: %s is not a valid source model for the Green''s function', ...
        upper(mfilename),src);
end
