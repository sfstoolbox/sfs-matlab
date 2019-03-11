function Ppwd = pwd_mono_circexp(Pm,Npw)
%PWD_MONO_CIRCEXP converts a circular basis expansion of a sound field to its
%two-dimensional plane wave decomposition
%
%   Usage: Ppwd = pwd_mono_circexp(Pm,[Npw])
%
%   Input parameters:
%       Pm      - circular basis expansion [N x (2*M+1)]
%       Npw     - number of equi-angular distributed plane waves, optional, 
%                 default: 2*M+1
%
%   Output parameters:
%       Ppwd    - plane wave decomposition [N x Npw]
%
%   See also: driving_function_mono_localwfs_sbl

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2019 SFS Toolbox Developers                             *
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
% https://sfs.readthedocs.io                            sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input parameters ==================================
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
isargmatrix(Pm);
M = (size(Pm,2)-1)./2;
if nargin == nargmin
  Npw = 2*M+1;
else
  isargpositivescalar(Npw);
end


%% ===== Computation ====================================================
% Implementation of
%                 ___
% _               \
% P(phipw, w) =   /__     P (w) i^m  e^(+i m phipw)
%               m=-M..M    m
% with
%
% phipw = n * 2*pi/Npw

Ppwd = inverse_cht(bsxfun(@times,Pm,1i.^(-M:M)),Npw);
