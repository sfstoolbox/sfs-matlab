function y = vector_norm(x,dim)
%VECTOR_NORM calculates the norm for the rows/columns of a matrix
%
%   Usage: y = vector_norm(x,dim)
%
%   Input parameters:
%       x   - matrix [n x m]
%       dim - dimension along the norm should be calculated
%
%   Output parameter:
%       y   - norm of the vectors [1 x m] or [n x 1]
%
%   VECTOR_NORM(x,dim) calculates the p-norm (with p=2) for the vectors
%   given within the matrix along the dimension dim.
%
%   See also: norm, vector_product

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
% NOTE: this is disabled due to performance issues in HRTF extrapolation, where
% this function is called multiple times.
%nargmin = 2;
%nargmax = 2;
%narginchk(nargmin,nargmax)
%isargmatrix(x)
%isargpositivescalar(dim)


%% ===== Computation =====================================================
y = sum(abs(x).^2,dim).^(1/2);
