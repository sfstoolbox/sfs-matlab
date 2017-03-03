function [x0,n0,w0] = rounded_box(t,ratio)
%ROUNDED_BOX computes position and direction vectors in 3D given a parameter 
%   indicating the position on the boundary of a box with rounded corners.
%
%   Usage: [x0,n0,w0] = rounded_box(t,ratio)
%
%   Input options:
%       t      - parameter indicating position on the boundary [1 x n]
%       ratio  - ratio between bending radius of the rounded corners and the 
%                half edge length of the rectangular bounding box (0,1) 
%
%   Output options:
%       x0     - positions [n x 3]
%       n0     - unit vector for boundaries normal vector [n x 3]
%       w0     - weights for integration [n x 1]
%
%   The box has rectangular bounding box in the horizontal plane with an edge 
%   length of 2 and its center at [0,0,0]. Choosing 0.0 for ratio leads to a
%   square, while 1.0 yields a circle.
%
%   The parameter t indicates the position on the boundary starting for t=0 
%   yielding x0=[1,0,0] and following the boundary in a counter-clockwise
%   manner. This function is periodic with respect to t with a period of 1.
%   Given two parameters t1 and t2=t1+delta, this function yields two positions
%   x01=rounded_box(t1, ratio) and x02=rounded_box(t2, ratio), whose distance 
%   ALONG THE BOUNDARY of the box is rounded_box(delta, ratio).
%
%   n0 is orthogonal to the boundary and is pointing inwards the box. However, 
%   for ratio=1 and t=1/8, 3/8, 5/8, ..., i.e a location in a corner of a
%   square, n0 has an angle of 45 degree to the adjacent edges of the box.
%
%   For a correct computation of the integration weights the elements have
%   to be sorted, either in ascending or descending order. It is further 
%   assumed, that the last element of t is the "left" neighbor of the first
%   element. The weights are computed on the distance of neighboring elements
%   ALONG THE BOUNDARY of the box.
%
% See also: secondary_source_positions

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

%% ===== Checking of input  parameters =======================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargpositivescalar(ratio);
isargvector(t);

if ratio > 1.0
  error('%s: ratio has to be between 0.0 and 1.0', upper(mfilename));
end

%% ===== Variables ===========================================================

N = length(t);

ratio_prime = (1.0 - ratio);

% the rounded box can be divided into four quarters, where each can be handled
% in the same manner and is rotated afterwards
rot_ind = mod( floor(t*4), 4 ) + 1;  % index for rotating
tseg = mod(t, 0.25)*4;  % map t to a single quarter

% length of a quarter of the rounded box
% (quarter circle + 2*linear segment);
l = pi/2*ratio + 2*ratio_prime;

circle = ratio_prime./l;  % value of t, where circular segment begins
circle_prime = 1.0 - circle;  % value t, where circular segment ends

%% ===== Positions & Orientations ============================================

x0 = zeros(N, 3);
n0 = zeros(N, 3);

%%%%% first linear segment %%%%%
seglin1 = ( tseg < circle );
% x0 = [1.0, t*l, 0]
x0(seglin1, 1) = 1.0;
x0(seglin1, 2) = tseg(seglin1)*l;
% n0 = [-1.0, 0, 0]
n0(seglin1, 1) = -1.0;

%%%%% second linear segment %%%%%
seglin2 = ( tseg > circle_prime );
% x0 = [(1.0-t)*l, 1.0, 0]
x0(seglin2 , 1) = ( 1.0 - tseg(seglin2) ) * l;
x0(seglin2 , 2) = 1.0;
% n0 = [0, -1.0, 0]
n0(seglin2, 2) = -1.0;

%%%%% circular segment %%%%%
segcirc = ~( seglin1 | seglin2 );
if ratio == 0
  phi = pi/4;
else
  phi = (tseg(segcirc) - circle)./(circle_prime - circle)*pi/2;
end
% x0 = circle_center + ratio*[cos(phi), sin(phi), 0]
x0(segcirc, 1) = ratio_prime + ratio*cos(phi);
x0(segcirc, 2) = ratio_prime + ratio*sin(phi);
% n0 = [-cos(phi), -sin(phi), 0]
n0(segcirc, 1) = -cos(phi);
n0(segcirc, 2) = -sin(phi);

% rotate quarters
rotx = [1, 0, 0; 0,-1, 0;-1, 0, 0; 0, 1, 0];
roty = [0, 1, 0; 1, 0, 0; 0,-1, 0;-1, 0, 0];

x0(:,1:2) = [sum( x0 .* rotx(rot_ind, :), 2), sum( x0 .* roty(rot_ind, :), 2)];
n0(:,1:2) = [sum( n0 .* rotx(rot_ind, :), 2), sum( n0 .* roty(rot_ind, :), 2)];

%% ===== Weights =============================================================

% Distance between adjacent samples. This takes the special geometry into
% account and is faster than using the more general secondary_source_distance()
% function.
dist = mod( t - t([N,1:N-1]), 1);  % distance to neighbor
dist = min( abs(dist), abs(1 - dist) );  % take periodicity into account
% Half of the distance to the "left" and half of the distance to the
% "right" neighbor
w0 = ( 0.5*dist + 0.5*dist([2:N,1]) )*4*l;
w0 = w0.';

end
