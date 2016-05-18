function z = convolution(x,y)
%CONVOLUTION convolves the signals x and y
%
%   Usage: z = convolution(x,y)
%
%   Input parameters:
%       x       - matrix/vector with signals as columns
%       y       - matrix/vector with signals as columns, note that only one of
%                 the signals can be a matrix
%
%   Output parameters:
%       z       - convolved signal
%
%   CONVOLUTION(x,y) convolves the signals given with x and y. One of the input
%   signals can be a matrix containing the signals as column vectors, the other
%   one has to be a column vector. The convolution is done in the frequency
%   domain and it is checked if we have only real signals to speed up the
%   calculation. The length of z is length(x)+length(y)-1.
%
%   See also: fft_real, ifft_real, fft, ifft

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Team                                   *
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


%% ===== Checking input parameters =======================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargmatrix(x,y);
% Check if only one of the inputs is a matrix
if all(size(x)>1) && all(size(y)>1)
    error('%s: Only one of the inputs can be multi-dimensional.', ...
        upper(mfilename));
end
% Ensure column vectors
if ~all(size(x)>1), x=column_vector(x); end
if ~all(size(y)>1), y=column_vector(y); end


%% ===== Computation =====================================================
% If one of the input signals is a matrix repmat the vector of the other signal
if all(size(x)>1)
    y = repmat(y,1,size(x,2));
elseif all(size(y)>1)
    x = repmat(x,1,size(y,2));
end
% Length of output signal
N = size(x,1)+size(y,1)-1;
% Convolve the signals in frequency domain
if isreal(x) && isreal(y)
    z = ifft_real(fft_real(fix_length(x,N)).*fft_real(fix_length(y,N)),N);
else
    z = ifft(fft(fix_length(x,N)).*fft(fix_length(y,N)));
end
