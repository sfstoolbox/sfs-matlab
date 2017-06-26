function [A,Phi] = inverse_cht(Am,Nphi)
%INVERSE_CHT computes the inverse circular harmonics transform (ICHT)
%
%   Usage: [A,Phi] = inverse_cht(Am,[Nphi])
%
%   Input parameters:
%       Am      - circular harmonics coefficients [N x (2*M+1)]
%       Nphi    - number of equi-angular distributed angles, for which the ICHT
%                 is computed, optional, default: 2*M+1
%
%   Output parameters:
%       A       - inverse circular harmonics transform [N x Nphi]
%       Phi     - corresponding angle of the ICHT [1 x Nphi]
%
%   See also: pwd_imp_circexp

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


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
isargmatrix(Am);
if nargin == nargmin
    Nphi = size(Am, 2);
else
    isargpositivescalar(Nphi);
end


%% ===== Computation ==================================================
M = (size(Am,2)-1)/2;
N = size(Am,1);

% Implementation of
%           ___
%           \
% A(phi) =  /__    A  e^(-j*m*n*2*pi/Nphi)
%         m=-M..M   m

% Spatial IFFT
A = zeros(N, Nphi);
% this handles cases where Nphi < M
for l=1:N
    A(l,:) = sum(buffer(Am(l,:),Nphi),2);
end
A = circshift(A,[0,-M]);  % m = 0, ..., M, ..., -M, ..., -1
A = ifft(A,[],2) * Nphi;  % IFFT includes factor 1/Nphi

% Axis corresponding to ICHT
if nargout>1
    Phi = 0:2*pi / Nphi:2*pi*(1-1/Nphi);
end
