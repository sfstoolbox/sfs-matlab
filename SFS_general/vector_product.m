function y = vector_product(x1,x2,dim)
%VECTOR_PRODUCT calculates the scalar product for two vectors within matrices
%
%   Usage: y = vector_product(x1,x2,dim)
%
%   Input parameters:
%       x1  - first matrix [n x m]
%       x2  - second matrix [n x m]
%       dim - dimension along the scalar product should be calculated
%
%   Output parameter:
%       y   - scalar product between the vectors [1 x m] or [n x 1]
%
%   VECTOR_PRODUCT(x1,x2,dim) calculates the scalar product between the vectors
%   given within the matrices along the dimension dim.
%
%   See also: vector_norm

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
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


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 3;
narginchk(nargmin,nargmax)
isargmatrix(x1,x2)
isargpositivescalar(dim)
if size(x1)~=size(x2)
    error('%s: the given matrices x1, x2 had to be of the same size.', ...
        upper(mfilename));
end


%% ===== Computation =====================================================
y = sum(x1.*x2,dim);
