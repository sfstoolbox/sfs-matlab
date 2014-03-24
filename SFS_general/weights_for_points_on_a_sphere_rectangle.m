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
 
%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
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
    % if the next theta is not equal to the previous one -> calculate new weight
    if theta(ii)~=theta(ii+1) 
        % calculate rectangle around the point with 
        % b*a = (r*cos(theta)d_phi) * (r*d_theta)
        a = (theta(ii)-theta(ii+1));
        b = (phi(ii+1))-(phi(ii));
        weights(ii) = a*b; % calculate area element
    else
        weights(ii) = weights(ii-1);
    end
end

% pole weight shouldn't be zero
if weights(1)==0
    weights(1) = weights(2);
end

w = 100*[weights,weights(end-1)]';
