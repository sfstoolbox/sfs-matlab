function Pm = circexp_convert_pwd(Ppwd, Mforce)
%CIRCEXP_CONVERT_PWD converts a 2D plane wave decomposition into regular
%circular expansion coefficients
%
%   Usage: Pm = circexp_convert_pwd(Ppwd, [Mforce])
%
%   Input parameters:
%       Ppwd          - 2D plane wave decomposition [Npw x Nfft]
%       Mforce        - maximum order of circular expansion, optional
%
%   Output parameters:
%       Pm            - regular circular expansion coefficients
%
%  Second dimension of Ppwd can be time or time-frequency axis.

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

%% ===== Checking of input parameters ==================================
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
isargmatrix(Ppwd);
Npw = size(Ppwd,1);
M = floor(Npw/2);
if nargin == nargmin
  Mforce = M;
  Mmin = M;
else
  isargpositivescalar(Mforce);
  Mmin = min(M, Mforce);
end

%% ===== Variables ======================================================
Nfft = size(Ppwd,2);

%% ===== Computation ====================================================
% Numerical Implementation of
%             m  1  / _
% P (w) = (-j)  --- | P(phi0, w) e^(-i m phipw) d phipw
%  m            2pi /
%
% using Discrete Fourier Transform
%                    ___
%             m  1   \        _ / 2pi    \   /   2pi     \
% P (w) = (-j)  ---  /__      P | ---, w | e^|-i --- n m |
%  m            Npw  n=0..Npw   \ Npw    /   \   Npw     /

Pm = fft(Ppwd,[],1)/Npw;  % fft does not include normalisation
% m = 0..N --> m = -[N/2] ... [N/2]
Pm = fftshift(Pm,1);
if mod(Npw, 2) == 0
  Pm = [Pm; Pm(1,:)];
end
Pm = bsxfun(@times, Pm, (-1j).^(-M:M).');

% force order
Ptmp = zeros(2*Mforce+1, Nfft);
Ptmp( Mforce-Mmin+1:Mforce+Mmin+1, :) = Pm( M-Mmin+1:M+Mmin+1, :);
Pm = Ptmp;
