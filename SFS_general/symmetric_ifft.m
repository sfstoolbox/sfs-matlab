function X = symmetric_ifft(varargin)
%SYMMETRIC_IFFT computes the inverse discrete Fourier transform of using a
%Fast Fourier Transform (FFT) algorithm.
%
%   Usage: X = symmetric_ifft(Y,[n,[dim]])
%
%   Input parameters:
%       Y       - input signal
%       n       - number of elements of X to use
%       dim     - dimension along the ifft should be calculated
%
%   Output parameters:
%       X       - inverse Fourier trafo of X
%
%   SYMMETRIC_IFFT(Y,n,dim) calculates the inverse discrete Fourier transform and
%   treats Y as if it were conjugate symmetric. Internally this calls
%   ifft(Y,n,dim,'symmetric') in the case of Matlab and real(ifft(Y,n,dim)) in the
%   case of Octave as there is no check for symmetry available.
%
%   See also: ifft, modal_weighting, interpolate_ir

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


% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 3;
narginchk(nargmin,nargmax);


%% ===== Calculation =====================================================
if isoctave
    X = real(ifft(varargin{:}));
else
    X = ifft(varargin{:},'symmetric');
end
