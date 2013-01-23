function [ir1,ir2,ir3,x0,desired_point] = findnearestneighbour_scalar(irs,phi,delta,r,number_of_neighbours,X0)
%FINDNEARESTNEIGHBOUR_SCALAR finds the 2 or 3 nearest neighbours of a given
%                            position (phi,theta,r)
%
%   Usage: [ir1,ir2,ir3,x0,desired_point] = findnearestneighbour_scalar
%                                           (irs,phi,delta,r,number_of_neighbours,X0)
%
%   Input parameters:
%       irs                  - data set containing the IRs
%       phi                  - azimuth angle of the desired point
%       delta                - elevation angle of the desired point
%       r                    - radius of the desired point
%       number_of_neighbours - determines how much neighbours are
%                              calculated. 2 for 2d case or 3 for 3d case 
%                              are possible inputs
%       X0                   - listener position 
%
%   Output parameters:
%       ir1,ir2,ir3      - IRs next to the position of the desired IR  
%       x0               - Matrix containing postions of ir1,ir2,ir3
%       desired_point    - position where the IR should be calculated
%   
%FINDNEARESTNEIGHBOUR_SCALAR(irs,phi,delta,r,number_of_neighbours,X0)
%   Finds the 2 or 3 nearest neighbours of a desired position by 
%   calculating the scalar product of the desired point and all
%   points where IRs available from the data set given by irs.
%   
%
%   see also: get_ir, findnearestneighbour_distance

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************    
%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 6;
error(nargchk(nargmin,nargmax,nargin));
if nargin == 5
    X0 = [0 0 0];
end

% turn on(=1) / off(=0) the debug mode for this file
debug = 1;

%% ===== Computation ====================================================
% calculate the x-,y-,z-coordinates for the desired angles(phi,delta)
% and radius r
[x,y,z] = sph2cart(phi,delta,r);
desired_point = [x,y,z];


% calculate the x-,y-,z-coordinates for all known IRs
[x,y,z] = sph2cart(irs.apparent_azimuth(1,:), ...
irs.apparent_elevation(1,:),irs.distance(1,:));
x0 = [x;y;z];
    
if ~any(irs.apparent_azimuth(1,1)-irs.apparent_azimuth(1,2:end)) && ~any(irs.apparent_elevation(1,1)-irs.apparent_elevation(1,2:end))

    [~,idx] = findnearestneighbour_distance(x0,desired_point,X0);
    x0 = [x0(:,idx),[0;0;0]];
    ir1 = [irs.left(:,idx(1)) irs.right(:,idx(1))]; 
    ir2 = [irs.left(:,idx(2)) irs.right(:,idx(2))];
    ir3 = [0 0 0];
    
else
% calculate the cos(angle) between desired point and all points of 
% the HRIR set 
    scalarproduct = (desired_point*x0)';
    for ii = 1:length(x0)
        norm_vec(ii,1) = norm(x0(:,ii));
    end
    cos_angle = scalarproduct./(norm_vec.*norm(desired_point));
% get the x nearest neighbours of the desired point
    [~,idx] = sort(cos_angle,'descend');
    idx = idx(1:number_of_neighbours);
    
    % store the previous found nearest neighbours, if you want to prove
    % them
    if debug
    save('NearestPoints','idx');
    end
    
% get the HRIR-positions corresponding to the maxima
    x0 = x0(:,idx);
    ir1 = [irs.left(:,idx(1)) irs.right(:,idx(1))];
    ir2 = [irs.left(:,idx(2)) irs.right(:,idx(2))];
    ir3 = [irs.left(:,idx(3)) irs.right(:,idx(3))]; 

end

% show the indices of the nearest neighbours. maybe desired for proving
% results
if debug 
fprintf('The indices of the found IRs with get_ir_new_version are: index1 = %d, index2 = %d\n, index3 = %d',idx(1),idx(2),idx(3));
fprintf('\n'); 
end
    
end
       



