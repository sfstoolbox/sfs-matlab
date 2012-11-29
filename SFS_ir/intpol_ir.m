function ir = intpol_ir(ir1,beta1,ir2,beta2,alpha)
%INTPOL_IR interpolates two given IRs for the given angle
%
%   Usage: ir = intpol_ir(ir1,beta1,ir2,beta2,alpha)
%
%   Input parameters:
%       ir1     - IR with lower angle
%       beta1   - angle of ir1 (rad)
%       ir2     - IR with bigger angle
%       beta2   - angle of ir2 (rad)
%       alpha   - angle of the desired IR (rad)
%
%   Output parameters:
%       ir      - IR for the given angle alpha (length(IR1),2)
%
%   INTPOL_IR(ir1,beta1,ir2,beta2,alpha) interpolates the two given IRs ir1 and
%   ir2 with their corresponding angles beta1 and beta2 for the given angle
%   alpha and returns an interpolated IR.
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
nargmin = 5;
nargmax = 5;
narginchk(nargmin,nargmax);


%% ===== Computation ====================================================

% Check if the given IR have the same angle
% in order to get the right interpolation dimension
if beta1==beta2
    error('%s: The angles of the two given IRs are the same!',upper(mfilename));
end

% Correct the given angles
beta1 = correct_azimuth(beta1);
beta2 = correct_azimuth(beta2);
alpha = correct_azimuth(alpha);

%if alpha==beta1 || alpha==beta2
%    error('%s: no interpolation needed for the given alpha value.',...
%        upper(mfilename));
%end
if length(ir1)~=length(ir2)
    error('%s: the given IRs have not the same length.',upper(mfilename));
end

% Linear interpolate the two given IRs
ir = ir1 + (ir2-ir1) ./ (beta2-beta1)*(alpha-beta1);
