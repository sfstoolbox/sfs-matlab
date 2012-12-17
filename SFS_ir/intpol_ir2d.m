function ir = intpol_ir2d(ir1,ir2,x0,desired_point)
%INTPOL_IR2D interpolates two given IRs for the given angle
%
%   Usage: ir = intpol_ir2d(ir1,ir2,x0,desired_point)
%
%   Input parameters:
%       ir1             - IR 1
%       ir2             - IR 2
%       x0              - Matrix containing the positions of ir1 and ir2
%       desired_point   - point at which the 3 IRs should be interpolated
%                         (in cartesian coordinates)
%
%   Output parameters:
%       ir      - IR for the desired position
%
%   INTPOL_IR2D(ir1,ir2,x0,desired_point) interpolates the two given IRs 
%   ir1 and ir2 with their corresponding angles beta1 and beta2 for the 
%   given angle alpha and returns an interpolated IR.
%
%   see also: get_ir, shorten_ir, read_irs

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
nargmin = 4;
nargmax = 4;
error(nargchk(nargmin,nargmax,nargin));

if length(ir1)~=length(ir2)
    error('%s: the given IRs have not the same length.',upper(mfilename));
end

%% ===== Computation ====================================================
%calculate the spherical coordinates for the 2 next neighbours
[phi1,theta1,r1] = cart2sph(x0(1,1),x0(2,1),x0(3,1));
[phi2,theta2,r2] = cart2sph(x0(2,1),x0(2,2),x0(2,3));
[alpha,beta,R]  = cart2sph(desired_point(1),desired_point(2),desired_point(3));


% Correct the given angles
phi1 = correct_azimuth(phi1);
phi2 = correct_azimuth(phi2);
theta1 = correct_elevation(theta1);
theta2 = correct_elevation(theta2);
alpha = correct_elevation(alpha);
beta = correct_elevation(beta);

% decide in which direction you have to interpolate
if phi1 == phi2 && theta1~=theta2
    
    phi1 = theta1;
    phi2 = theta2;
    alpha = beta;

elseif phi1 == phi2 && theta1 == theta2
    
    phi1 = r1;
    phi2 = r2;
    alpha = R;
    
end

% Linear interpolate the two given IRs
ir = ir1 + (ir2-ir1) ./ (phi2-phi1)*(alpha-phi1);