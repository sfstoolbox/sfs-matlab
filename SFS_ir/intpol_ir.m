function ir = intpol_ir(varargin)
%INTPOL_IR interpolates three given IRs for the given angle
%
%   Usage: ir = intpol_ir(ir1,ir2,[ir3],x0,xs)
%
%   Input parameters:
%       ir1     - IR 1
%       ir2     - IR 2
%       ir3     - IR 3 (optional)
%       x0      - matrix containing positions of single IRs / rad
%       xs      - desired position after interpolation / rad
%
%   Output parameters:
%       ir      - IR for the given position
%
%   INTPOL_IR(ir1,phi1,theta1,ir2,phi2,theta2,ir3,phi3,theta3,alpha,beta)
%   interpolates the three given IRs ir1,ir2 and ir3 with their corresponding 
%   angles (phi1,theta1),(phi2,theta2) and (phi3,theta3) for the given
%   angles (alpha,beta) and returns an interpolated IR.
%   Note that the given parameter are not checked if they have all the right
%   dimensions in order to save computational time, because this function could
%   be called quiet often.
%
%   see also: get_ir, shorten_ir, read_irs

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
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);

%% ===== Computation ====================================================
% --- 1D interpolation ---
if nargin==4
    ir1 = varargin{1};
    ir2 = varargin{2};
    x0 = varargin{3};
    xs = varargin{4};
    % linear interpolation
    ir = ir1 + (ir2-ir1) * norm(xs-x0(:,1))/norm(x0(:,2)-x0(:,1));
else
    ir1 = varargin{1};
    ir2 = varargin{2};
    ir3 = varargin{3};
    x0 = varargin{4};
    xs = varargin{5};
    % linear interpolation
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
    ir = w(1)*ir1 + w(2)*ir2 + w(3)*ir3;
end
