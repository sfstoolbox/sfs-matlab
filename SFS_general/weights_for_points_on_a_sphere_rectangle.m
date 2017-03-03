function [w] = weights_for_points_on_a_sphere_rectangle(phi,theta)
%WEIGHTS_FOR_POINTS_ON_A_SPHERE_RECTANGLE returns the weights for a given
% set of points on a sphere.
%
%   Usage: [w] = weights_for_points_on_a_sphere_rectangle(phi,theta)
%
%   Input parameters:
%       phi           - azimuth angles of given points (row vector) / rad
%       theta         - elevation angles of given points (row vector) / rad
%
%   Output parameters:
%       w             - weights for given set of points on the sphere
%
%   WEIGHTS_FOR_POINTS_ON_A_SPHERE_RECTANGLE(phi,theta) returns the weights
%   for a given set of points on a sphere. The weights will be achieved by
%   calculating a rectangle around every point given by its azimuth angle
%   (phi), elevation (theta) angle.
%
%   See also: get_spherical_grid

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
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);


%% ===== Computation ====================================================
weights = zeros(1,length(theta)-1);
% Start with first element (supposed to be a pole)
% Pole d_phi to the next point would be inf. Therefore the 2nd nearest d_phi is
% used at this time
a = (theta(1)-theta(2));
b = (phi(4))-(phi(3));
weights(1) = a*b;
for ii=2:length(theta)-1
    % If the next theta is not equal to the previous one -> calculate new weight
    if theta(ii)~=theta(ii+1) 
        % Calculate rectangle around the point with 
        % b*a = (r*cos(theta)d_phi) * (r*d_theta)
        a = (theta(ii)-theta(ii+1));
        b = (phi(ii+1))-(phi(ii));
        weights(ii) = a*b; % calculate area element
    else
        weights(ii) = weights(ii-1);
    end
end

% Pole weight shouldn't be zero
if weights(1)==0
    weights(1) = weights(2);
end

w = 100*[weights,weights(end-1)]';
