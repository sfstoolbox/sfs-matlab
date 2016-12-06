function B = interpolation(A,X,x)
%INTERPOLATION interpolates the given data A(X) at the point x
%
%   Usage: B = interpolation(A,X,x)
%
%   Input parameters:
%       A       - matrix containing data as rows in the form [N M], where
%                     N ... samples of data A
%                     M ... number of points X (2 or 3)
%       X       - matrix containing positions b as columns, at which A(X) is
%                 given [D M]
%                     D ... dimension of space (1 or 2)
%                     M ... number of points X (2 or 3)
%       x       - desired point at which A should be interpolated [D 1]
%                     D ... dimension of space (1 or 2)
%
%   Output parameters:
%       B       - interpolated data at point x [N 1]
%
%   INTERPOLATION(A,X,x) linearly interpolates the data A given at two or three
%   points X at the desired position x.
%
%   See also: interpolate_ir, get_ir

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input parameters ===================================
% Disabled for time constraints
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
    B = (w(1)*A(:,1) + w(2)*A(:,2) + w(3)*A(:,3)) / sum(w);
else
    error('%s: size(A,2) has to be 2 or 3.',upper(mfilename));
end
