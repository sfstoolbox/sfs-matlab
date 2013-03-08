function ir = intpol_ir(desired_point,x0_1,ir1,x0_2,ir2,x0_3,ir3)
%INTPOL_IR interpolates three given IRs for the given angle
%
%   Usage: ir = intpol_ir(ir1,phi1,theta1,ir2,phi2,theta2,ir3,phi3,theta3,alpha,beta)
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
%   INTPOL_IR(ir1,phi1,theta1,ir2,phi2,theta2,ir3,phi3,theta3,alpha,beta)
%   interpolates the three given IRs ir1,ir2 and ir3 with their corresponding 
%   angles (phi1,theta1),(phi2,theta2) and (phi3,theta3) for the given
%   angles (alpha,beta) and returns an interpolated IR.
%
%   see also: get_ir, shorten_ir, read_irs

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 7;
error(nargchk(nargmin,nargmax,nargin));

if nargin == 5

    if length(ir1)~=length(ir2)
    error('%s: the given IRs have not the same length.',upper(mfilename));
    end
    
    L = [x0_1,x0_2];
    p = desired_point';
    g = p\L^-1;
    % calculate desired ir with linear combination of ir1,ir2 
    ir = g(1,1)*ir1 + g(1,2)*ir2;
    

%% ===== Computation ====================================================
else
    % check if the length of the found IRs are the same
    if length(ir1)~=length(ir2) || length(ir2)~=length(ir3) || length(ir1)~=length(ir3)
        error('%s: the given IRs have not the same length.',upper(mfilename));
    end

    % Solve linear equation system to get the desired weight factors g(n)
    L = [x0_1,x0_2,x0_3];
    p = desired_point';
    g = p\L^-1;
    % calculate desired ir with linear combination of ir1,ir2 and ir3
    ir = g(1,1)*ir1 + g(1,2)*ir2 + g(1,3)*ir3;
    
end
