function ir_new = intpol_ir(ir,x0,xs)
%INTPOL_IR interpolates three given IRs for the given angle
%
%   Usage: ir = intpol_ir(ir,x0,xs)
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
%
%   INTPOL_IR(ir,x0,xs)
%   interpolates the two to three given impulse responses from ir with their
%   corresponding angles x0 for the given angles xs and returns an interpolated
%   impulse response.
%   Note that the given parameter are not checked if they have all the correct
%   dimensions in order to save computational time, because this function could
%   be called quiet often.
%
%   See also: get_ir, shorten_ir, SOFAload

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
nargmin = 3;
nargmax = 3;
narginchk(nargmin,nargmax);


%% ===== Computation ====================================================
% --- 1D interpolation ---
if size(ir,1)==2
    % Linear interpolation
    ir_new = ir(1,:,:) + (ir(2,:,:)-ir(1,:,:)) * ...
         norm(xs-x0(:,1)) / norm(x0(:,2)-x0(:,1));
% --- 2D interpolation ---
elseif size(ir,1)==3
    % Linear interpolation (compare Vector Based Amplitude Panning)
    %
    %           x0(:,ii) xs
    % w(ii) = --------------
    %         |x0(:,ii)||xs|
    %
    w = vector_product(x0,repmat(xs',[1 3]),1) / ...
        ( vector_norm(x0,1)./norm(xs));
    % The interpolation with 3 points hasn't been checked yet, hence we are
    % including a checking of the w parameters
    if any(w<0)
        error('%s: one of your interpolation weights is <0.',upper(mfilename));
    end
    % calculate desired ir with linear combination of ir1,ir2 and ir3
    ir_new = w(1)*ir(1,:,:) + w(2)*ir(2,:,:) + w(3)*ir(3,:,:);
else
    error('%s: size(ir,1) has to be 2 or 3.',upper(mfilename));
end
