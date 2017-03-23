function G = greens_function_mono(x,y,z,xs,src,f,conf)
%GREENS_FUNCTION_MONO returns a Green's function in the frequency domain
%
%   Usage: G = greens_function_mono(x,y,z,xs,src,f,conf)
%
%   Input options:
%       x,y,z   - x,y,z points for which the Green's function should be
%                 calculated / m
%       xs      - position of the source
%       src     - source model of the Green's function. Valid models are:
%                   'ps'  - point source
%                   'ls'  - line source
%                   'pw'  - plane wave
%                   'dps' - dipole point source
%       f       - frequency of the source / Hz
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       G       - Green's function evaluated at the points x,y,z
%
%   GREENS_FUNCTION_MONO(x,y,z,xs,src,f,conf) calculates the Green's function
%   for the given source model located at xs for the given points x,y and the
%   frequency f.
%
%   See also: sound_field_mono

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
% Disabled checking for performance reasons


%% ===== Configuration ==================================================
c = conf.c;
phase = conf.phase;


%% ===== Computation =====================================================
% Frequency
omega = 2*pi*f;
% Calculate Green's function for the given source model
if strcmp('ps',src)
    % Source model for a point source: 3D Green's function.
    %
    %              1  e^(-i w/c |x-xs|)
    % G(x-xs,w) = --- -----------------
    %             4pi      |x-xs|
    %
    % See http://sfstoolbox.org/#equation-S.ps
    %
    G = 1/(4*pi) * exp(-1i*omega/c .* sqrt((x-xs(1)).^2+(y-xs(2)).^2+(z-xs(3)).^2)) ./ ...
            sqrt((x-xs(1)).^2+(y-xs(2)).^2+(z-xs(3)).^2);

elseif strcmp('dps',src)
    % Source model for a dipole point source: derivative of 3D Green's function.
    %
    %  d                1   / iw       1    \   (x-xs) ns
    % ---- G(x-xs,w) = --- | ----- + ------- | ----------- e^(-i w/c |x-xs|)
    % d ns             4pi  \  c     |x-xs| /   |x-xs|^2
    %
    % See http://sfstoolbox.org/#equation-S.dps
    %
    % r = |x-xs|
    r = sqrt((x-xs(1)).^2+(y-xs(2)).^2+(z-xs(3)).^2);
    % scalar = (x-xs) nxs
    scalar = xs(4).*(x-xs(1)) + xs(5).*(y-xs(2))  + xs(6).*(z-xs(3));
    %
    G = 1/(4*pi) .* (1i*omega/c + 1./r) .* scalar./r.^2 .* exp(-1i*omega/c.*r);

elseif strcmp('ls',src)
    % Source model for a line source: 2D Green's function.
    %
    %                i   (2) / w        \
    % G(x-xs,w) =  - -  H0  |  - |x-xs|  |
    %                4       \ c        /
    %
    % See http://sfstoolbox.org/#equation-S.ls
    %
    G = -1i/4 * besselh(0,2,omega/c* ...
        sqrt( (x-xs(1)).^2 + (y-xs(2)).^2 + (z-xs(3)).^2 ));

elseif strcmp('pw',src)
    % Source model for a plane wave:
    %
    % G(x,w) = e^(-i w/c n x)
    %
    % See: http://sfstoolbox.org/#equation-S.pw
    %
    % Direction of plane wave
    nxs = xs(:,1:3) / norm(xs(:,1:3));
    % Calculate sound field
    G = exp(-1i*omega/c.*(nxs(1).*x+nxs(2).*y+nxs(3).*z));
else
    error('%s: %s is not a valid source model for the Green''s function', ...
        upper(mfilename),src);
end

% Add phase to be able to simulate different time steps
G = G .* exp(-1i*phase);
