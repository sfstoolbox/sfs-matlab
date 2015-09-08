function [x0, n0] = smooth_box(t, ratio)
%SMOOTH_BOX computes position and direction vectors in 3D given a parameter 
%   indicating the position on the boundary of a box with smooth corners.
%
%   Usage: [x0, n0] = smooth_box(t, ratio)
%
%   Input options:
%       t      - parameter indicating position on the boundary [1 x n]
%       ratio  - ratio between bending radius of smooth corners and the 
%                half edge length of the rectangular bounding box (0,1) 
%
%   Output options:
%       x0     - positions [n x 3]
%       n0     - unit vector for boundaries normal vector [n x 3]
%
%   The box has rectangular bounding box in the horizontal plane with an edge 
%   length of 2 and its center at [0,0,0]. Choosing 0.0 for ratio leads to a
%   square with sharp corners,  while 1.0 yields a circle. 
%   The parameter t indicates the position on the  boundary starting for t=0 
%   yielding x0=[1,0,0] and following the boundary in a counter-clockwise
%   manner. This function is period with respect to t with a period of 1.
%   Given two parameters t1 and t2=t1+delta, this function yields two positions
%   x01=smooth_box(t1, ratio) and x02=smooth_box(t2, ratio), whose distance 
%   ALONG THE BOUNDARY of the box is smooth_box(delta, ratio).
%   n0 is orthogonal to the boundary. However, for ratio=1 and t=1/8, 3/8, 
%   5/8, ..., i.e a location in a sharp corner, n0 has an angle of 45 degree
%   to the adjacent edges of the box.
%
% See also: secondary_source_positions

%*****************************************************************************
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
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

%% ===== Computation =========================================================

ratio_prime = (1.0 - ratio);

% the smooth box can be divided into four quarters, where each can be handled
% in the same manner and is rotated afterwards
rot_ind = mod( floor(t*4), 4 ) + 1;  % index for rotating
t = mod(t, 0.25)*4;  % map t to a single quarter

% length of a quarter of the smoothed box
% (quarter circle + 2*linear segment);
l = pi/2*ratio + 2*ratio_prime;

circle = ratio_prime./l;  % value of t, where circular segment begins
circle_prime = 1.0 - circle;  % value t, where circular segment ends

x0 = zeros(length(t), 3);
n0 = zeros(length(t), 3);

%%%%% first linear segment %%%%%
seglin1 = ( t < circle );
% x0 = [1.0, t*l, 0]
x0(seglin1, 1) = 1.0;
x0(seglin1, 2) = t(seglin1)*l;
% n0 = [-1.0, 0, 0]
n0(seglin1, 1) = -1.0;

%%%%% second linear segment %%%%%
seglin2 = ( t > circle_prime );
% x0 = [(1.0-t)*l, 1.0, 0]
x0(seglin2 , 1) = ( 1.0 - t(seglin2) ) * l;
x0(seglin2 , 2) = 1.0;
% n0 = [0, -1.0, 0]
n0(seglin2, 2) = -1.0;

%%%%% circular segment %%%%%
segcirc = ~( seglin1 | seglin2 );
if ratio == 0
  phi = pi/4;
else
  phi = (t(segcirc) - circle)./(circle_prime - circle)*pi/2;
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

end
