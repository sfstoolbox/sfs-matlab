function [ir_new,x0_new] = interpolate_ir(ir,x0,xs,conf)
%INTERPOLATE_IR interpolates three given IRs for the given angle
%
%   Usage: ir = interpolate_ir(ir,x0,xs)
%
%   Input parameters:
%       ir      - matrix containing impulse responses in the form [M C N], where
%                     M ... Number of measurements (2<=M<=3)
%                     C ... Number of channels
%                     N ... Number of samples
%       x0      - matrix containing positions of single impulse
%                 responses [2 M] / (rad, rad)
%       xs      - desired position after interpolation [2 1] / (rad, rad)
%
%   Output parameters:
%       ir      - impulse response for the given position [1 C N]
%       x0      - position corresponding to the returned impulse response
%
%   INTERPOLATE_IR(ir,x0,xs)
%   interpolates the two to three given impulse responses from ir with their
%   corresponding angles x0 for the given angles xs and returns an interpolated
%   impulse response.
%   Note that the given parameter are not checked if they have all the correct
%   dimensions in order to save computational time, because this function could
%   be called quiet often.
%
%   See also: get_ir, interpolation

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


%% ===== Checking of input parameters ===================================
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);


%% ===== Configuration ==================================================
useinterpolation = conf.ir.useinterpolation;
% Precission of the wanted angle. If an impulse response within the given
% precission could be found no interpolation is applied.
prec = 0.001; % ~ 0.05 deg


%% ===== Computation ====================================================
ir_new = ir(1,:,:);
x0_new = xs;
% Check if we have found directly the desired point or have to interpolate
% bewteen different impulse responses
if norm(x0(:,1)-xs)<prec || ~useinterpolation
    % Return the first nearest neighbour
    x0_new = x0(:,1);
    return;
else
    % === IR interpolation ===
    % Check if we have to interpolate in one or two dimensions
    if norm(x0(1,1)-x0(1,2))<prec || norm(x0(2,1)-x0(2,2))<prec
        % --- 1D interpolation ---
        warning('SFS:irs_intpol',...
            ['doing 1D IR interpolation between (%.1f,%.1f) deg ', ...
             'and (%.1f,%.1f) deg.'], ...
            deg(x0(1,1)), deg(x0(2,1)), ...
            deg(x0(1,2)), deg(x0(2,2)));
        ir_new(1,1,:) = interpolation(squeeze(ir(1:2,1,:))',x0(1:2,:),xs);
        ir_new(1,2,:) = interpolation(squeeze(ir(1:2,2,:))',x0(1:2,:),xs);
    else
        % --- 2D interpolation ---
        warning('SFS:irs_intpol3D',...
            ['doing 2D IR interpolation between (%.1f,%.1f) deg, ', ...
             '(%.1f,%.1f) deg and (%.1f,%.1f) deg.'], ...
            deg(x0(1,1)), deg(x0(2,1)), ...
            deg(x0(1,2)), deg(x0(2,2)), ...
            deg(x0(1,3)), deg(x0(2,3)));
        ir_new(1,1,:) = interpolation(squeeze(ir(1:3,1,:))',x0(1:3,:),xs);
        ir_new(1,2,:) = interpolation(squeeze(ir(1:3,2,:))',x0(1:3,:),xs);
    end
end
