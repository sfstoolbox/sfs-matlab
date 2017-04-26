function ppwd = pwd_imp_circexp(pm,Npw)
%PWD_IMP_CIRCEXP converts a circular basis expansion of a sound field to its
%two-dimensional plane wave decomposition
%
%   Usage: ppwd = pwd_imp_circexp(pm,[Npw])
%
%   Input parameters:
%       pm      - circular basis expansion [N x (M+1)]
%       Npw     - number of equi-angular distributed plane waves, optional, 
%                 default: 2*M+1
%
%   Output parameters:
%       ppwd    - plane wave decomposition [N x Npw]
%
%   See also: driving_function_imp_localwfs_sbl_ps,
%   driving_function_imp_localwfs_sbl_pw

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
isargmatrix(pm);
M = size(pm,2)-1;
if nargin == nargmin
  Npw = 2*M+1;
else
  isargpositivescalar(Npw);
end


%% ===== Computation ====================================================
% Implementation of
%                 ___
% _               \
% p(phipw, t) =   /__     p (t) j^m  e^(-j m phipw)
%               m=-M..M    m
% with
%
% phipw = n * 2*pi/Npw

pm = [conj(pm(:,end:-1:2)), pm];  % append coefficients for negative m
ppwd = inverse_cht(bsxfun(@times,pm,1j.^(-M:M)),Npw);
