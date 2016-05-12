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
%       t       - time / samples
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       g       - Green's function evaluated at the points x,y,z
%       t       - Correspondiong time values with integrated time
%                 shift / samples
%
%   GREENS_FUNCTION_IMP(x,y,z,xs,src,t,conf) calculates the Green's function for
%   the given source model located at xs for the given points x,y,z. Note, that
%   the delta function for the time t is returned as an extra argument. If you
%   want the value of the Green's function only to this specific time you should
%   have a look at sound_field_imp() and apply the folowing command:
%   [p,x,y,z] = sound_field_imp(X,Y,Z,[xs 0 -1 0],src,1,t,conf);
%
%   References:
%       H. Wierstorf, J. Ahrens, F. Winter, F. Schultz, S. Spors (2015) -
%       "Theory of Sound Field Synthesis"
%
%   See also: greens_function_mono, sound_field_imp

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Team                                   *
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
fs = conf.fs;


%% ===== Computation =====================================================
% calculate Green's function for the given source model
if strcmp('ps',src)
    % Source model for a point source: 3D Green's function.
    %
    %                  1
    % g(x-xs,t) = ---------- delta(t - |x-xs|/c)
    %             4pi |x-xs|
    %
    % see Wierstorf et al. (2015), eq.(#s:ps)
    %
    r = sqrt((x-xs(1)).^2+(y-xs(2)).^2+(z-xs(3)).^2);
    g = 1./(4*pi.*r);
    t = (r/c)*fs-t;

elseif strcmp('dps',src)
    % Source model for a dipole point source: derivative of 3D Green's function.
    %
    %                 1   / -1 / iw \      1    \  (x-xs)ns
    % g(x-xs,ns,t) = --- | F  | ---- | + ------  | --------- delta(t - |x-xs|/c)
    %                4pi  \    \ c  /    |x-xs| /  |x-xs|^2
    %
    % see Wierstorf et al. (2015), eq.(#s:dps)
    %
    to_be_implemented(mfilename);

elseif strcmp('ls',src)
    % Source model for a line source: 2D Green's function.
    %                          ___
    %              -1/ c\     | 1       1
    % g(x-xs,t) = F |--  |  - |---  --_-_-_- delta(t - |x-xs|/c)
    %                \iw/    \|8pi  \||x-xs|
    %
    % see Wierstorf et al. (2015), eq.(#s:ls)
    % Note, that the filter F^-1 is not implemented!!!!
    %
    r = sqrt((x-xs(1)).^2+(y-xs(2)).^2+(z-xs(3)).^2);
    g = 1./sqrt(r) * sqrt(1/(8*pi));
    t = (r/c)*fs-t;

elseif strcmp('pw',src)
    % Source model for a plane wave:
    %
    % g(x,t) = delta(t - nx/c)
    %
    % see Wierstorf et al. (2015), eq.(#s:pw)
    %
    % direction of plane wave
    nxs = xs / norm(xs);
    %
    % The following code enables us to replace this two for-loops
    % for ii = 1:size(x,1)
    %     for jj = 1:size(x,2)
    %         t(ii,jj) = nxs*[x(ii,jj) y(ii,jj) z(ii,jj)]'./c;
    %     end
    % end
    %
    % Get a matrix in the form of
    % 1 1 1 0 0 0 0 0 0
    % 0 0 0 1 1 1 0 0 0
    % 0 0 0 0 0 0 1 1 1
    E = eye(3*size(x,1));
    E = E(1:3:end,:)+E(2:3:end,:)+E(3:3:end,:);
    % Multiply this matrix with the plane wave direction
    N = repmat(nxs,size(x,1)) .* E;
    % Interlace x,y,z into one matrix
    % x11 x12 ... x1m
    % y11 y12 ... y1m
    % z11 z12 ... z1m
    % .   .       .
    % .   .       .
    % xn1 xn2 ... xnm
    % yn1 yn2 ... ynm
    % zn1 zn2 ... znm
    XYZ = zeros(3*size(x,1),size(x,2));
    XYZ(1:3:end,:) = x;
    XYZ(2:3:end,:) = y;
    XYZ(3:3:end,:) = z;
    %
    g = 1;
    t = N*XYZ./c*fs-t;

else
    error('%s: %s is not a valid source model for the Green''s function', ...
        upper(mfilename),src);
end
