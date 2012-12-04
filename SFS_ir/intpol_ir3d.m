function ir = intpol_ir3d(ir1,phi1,theta1,ir2,phi2,theta2,ir3,phi3,theta3,alpha,beta)
%INTPOL_IR3D interpolates three given IRs for the given angle
%
%   Usage: ir = intpol_ir3d(ir1,phi1,theta1,ir2,phi2,theta2,ir3,phi3,theta3,alpha,beta)
%
%   Input parameters:
%       ir1     - IR 1
%       phi1    - azimuth angle of ir1 (rad)
%       theta1  - elevation angle of ir1 (rad)
%       ir2     - IR 2
%       phi2    - azimuth angle of ir2 (rad)
%       theta2  - elevation angle of ir2 (rad)
%       ir3     - IR 3
%       phi3    - azimuth angle of ir3 (rad)
%       theta3  - elevation angle of ir3 (rad)
%       alpha   - azimuth angle of the desired IR (rad)
%       beta    - elevation angle of the desired IR (rad)
%
%   Output parameters:
%       ir      - IR for the given angles alpha,beta (length(IR1),2)
%
%   INTPOL_IR3d(ir1,phi1,theta1,ir2,phi2,theta2,ir3,phi3,theta3,alpha,beta)
%   interpolates the three given IRs ir1,ir2 and ir3 with their corresponding 
%   angles (phi1,theta1),(phi2,theta2) and (phi3,theta3) for the given
%   angles (alpha,beta) and returns an interpolated IR.
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
nargmin = 11;
nargmax = 11;
error(nargchk(nargmin,nargmax,nargin));


%% ===== Computation ====================================================


% note: phi_min < alpha < phi_max AND theta_min < beta < theta_max
% otherwise you're not interpolating between them
phi = [phi1,phi2,phi3];
theta = [theta1,theta2,theta3];
if alpha < min(phi) || alpha > max(phi)
    error('%s: The azimuth angle of the desired IR has to be between phi1,phi2 and phi3!',upper(mfilename));
end
if beta < min(theta) || beta > max(theta)
    error('%s: The elevation angle of the desired IR has to be between theta1,theta2 and theta3!',upper(mfilename));
end
% Correct the given angles
phi1 = correct_azimuth(phi1);
phi2 = correct_azimuth(phi2);
alpha = correct_azimuth(alpha);

theta1 = correct_elevation(theta1);
theta2 = correct_elevation(theta2);
beta = correct_elevation(beta);

if length(ir1)~=length(ir2) || length(ir2)~=length(ir3) || length(ir1)~=length(ir3)
    error('%s: the given IRs have not the same length.',upper(mfilename));
end

% Linear interpolate the three given IRs
% calculate scaling parameters
m = (alpha - phi1)./(phi2 - phi1);
n = (beta - theta1)./(theta2 - theta1);

% Check if the given IR have the same angle
% in order to get the right interpolation dimension
if phi1==phi2 && phi2==phi3
    m=0;
    warning('%s: The azimuth angles of the three given IRs are the same!',upper(mfilename));
end

if theta1==theta2 && theta2==theta3
    n = 0;
    warning('%s: The elevation angles of the three given IRs are the same!',upper(mfilename));
end

% calculate desired ir
ir = ir1 + m.*(ir2-ir1) + n.*(ir3-ir1);
