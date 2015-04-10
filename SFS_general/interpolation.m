function B = interpolation(A,X,x)
%INTPOLATION interpolates the given data A(X) at the point x
%
%   Usage: B = interpolation(A,X,x)
%
%   Input parameters:
%       A       - matrix containing data as rows in the form [N M], where
%                     M ... number of points X (2 or 3)
%                     N ... samples of data A
%       X       - matrix containing positions b as columns, at which A(X) is
%                 given [D M]
%                     M ... number of points X (2 or 3)
%                     D ... dimension of space (1 or 2)
%       x       - desired point at which A should be interpolated [D 1]
%
%   Output parameters:
%       B       - interpolated data at point x [N 1]
%
%   INTERPOLATION(A,X,x) interpolates the data A given at two or three points X
%   at the desired position x.
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
% Disabled for time constrains
%nargmin = 3;
%nargmax = 3;
%narginchk(nargmin,nargmax);


%% ===== Computation ====================================================
% --- 1D interpolation ---
if size(A,2)==2
    % Linear interpolation
    B = A(:,1) + (A(:,2)-A(:,1)) * ...
         norm(x-X(:,1)) / norm(X(:,2)-X(:,1));
% --- 2D interpolation ---
elseif size(A,2)==3
    % Interpolation between three points (compare Vector Based Amplitude Panning)
    %
    %           X(:,ii) x
    % w(ii) = ------------
    %         |X(:,ii)||x|
    %
    % FIXME: this is not working if |X(:,ii)| or |x| = 0!
    w = vector_product(X,repmat(x,[1 3]),1) ./ ...
        (vector_norm(X,1)./norm(x));
    % The interpolation with 3 points hasn't been checked yet, hence we are
    % including a checking of the w parameters
    if any(w<0)
        error('%s: one of your interpolation weights is <0.',upper(mfilename));
    end
    % Calculate desired B with linear combination w(ii)
    B = w(1)*A(:,1) + w(2)*A(:,2) + w(3)*A(:,3);
else
    error('%s: size(A,2) has to be 2 or 3.',upper(mfilename));
end
