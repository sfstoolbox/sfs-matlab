function z = convolution(x,y)
%CONVOLUTION convolves the signals x and y
%
%   Usage: z = convolution(x,y)
%
%   Input parameters:
%       x       - matrix/vector with signals as columns
%       y       - matrix/vector with signals as columns
%
%   Output parameters:
%       z       - convolved signal
%
%   CONVOLUTION(x,y) convolves the signals given with x and y. If both input
%   signals are matrices, they must contain the same number of signals (columns).
%   The convolution is done in the frequency domain and it is checked if we have
%   only real signals to speed up the calculation. The length of z is
%   length(x)+length(y)-1.
%
%   See also: fft, ifft

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


%% ===== Checking input parameters =======================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargmatrix(x,y);
% Ensure column vectors
if ~all(size(x)>1), x=x(:); end
if ~all(size(y)>1), y=y(:); end
% If the inputs are two matrices, check if the number of signals is the same
if size(x,2)>1 && size(y,2)>1 && size(x,2)~=size(y,2)
  error(['%s: Two input matrices must have the same number of signals '...
      '(columns).'],upper(mfilename));
end

%% ===== Computation =====================================================
% Length of output signal
N = size(x,1)+size(y,1)-1;
% FFT + multiplication
Z = bsxfun(@times,fft(x,N,1),fft(y,N,1));  % automatically adjusts sizes
% IFFT
if isreal(x) && isreal(y)
    z = real(ifft(Z,[],1));
else
    z = ifft(Z,[],1);
end
