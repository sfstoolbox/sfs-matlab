function [tau, R, Delta] = retarded_time(x, y, z, t, vs, conf)
% RETARDED_TIME computes the retarded time for a uniformly moving point source
%
%   Usage: [tau, R, Delta] = retarded_time(x, y, z, t, vs, conf)
%
%   Input parameters:
%       x           - x coordinate position of observer / m [nx1] or [1x1]
%       y           - y coordinate position of observer / m [nx1] or [1x1]
%       z           - z coordinate position of observer / m [nx1] or [1x1]
%       t           - absolute time / s [nxm] or [1xm]
%       vs          - velocity vector of sound source / (m/s) [1x3] or [nx3]
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       tau         - retarded time between sound source and observer / s [nxm]
%       R           - retarded distance (R = tau*c) / m [nxm]
%       Delta       - auxilary distance / m [nxm]
%
%   RETARDED_TIMES(x, y, z, t, vs, conf) the retarded time for a point source
%   uniformly moving along a line. Hereby, vs defines its direction
%   and its velocity. The position of the sound source at t=0 is the coordinate
%   origin.
%
%   See also:

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
% nargmin = 6;
% nargmax = 6;
% narginchk(nargmin,nargmax);
% isargmatrix(x,y,z,t,vs);
% isargstruct(conf);


%% ===== Configuration ==================================================
c = conf.c;


%% ===== Computation ====================================================

v = vector_norm(vs,2);  % velocity of sound source [n x 1] or [1 x 1]
nvs = bsxfun(@rdivide,vs,v);  % direction of movement [n x 3] or [1 x 3]
M = v/c;  % Mach number [n x 1] or [1 x 1]

% vector between observer and point source d(t) = x - vs*t
% [n x m] or [1 x m]
dx = bsxfun(@minus, x, bsxfun(@times,t, vs(:,1)));
dy = bsxfun(@minus, y, bsxfun(@times,t, vs(:,2)));
dz = bsxfun(@minus, z, bsxfun(@times,t, vs(:,3)));
% instant distance between sound source and observer position
r = sqrt( dx.^2 + dy.^2 + dz.^2 );  % [n x m]
% component of d(t) parallel to the direction of movement: d(t)^T * nvs
% [n x m]
xpara = nvs(:,1).*dx + nvs(:,2).*dy + nvs(:,3).*dz;
% auxiliary distance  
Delta = sqrt( M.^2.*xpara.^2 + (1-M.^2).*r.^2 );  % [n x 1]
% retarded distance
R = (M.*xpara + Delta)./(1-M.^2);
% retarded time
tau = R./c;

end
