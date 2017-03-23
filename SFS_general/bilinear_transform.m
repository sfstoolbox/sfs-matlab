function [b,a] = bilinear_transform(sos,conf)
%BILINEAR_TRANSFORM returns the second-order section as filter coefficients
%
%   Usage: [b,a] = bilinear_transform(sos,conf)
%
%   Input parameters:
%       sos     - second-order section representation
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       b,a     - filter coefficients as cells
%
%   BILINEAR_TRANSFORM(sos,conf) transforms the second-order section
%   representation of a system to time discrete filter coefficients.
%
%   See also: driving_function_imp_nfchoa, zp2sos

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
isargmatrix(sos);
isargstruct(conf);


%% ===== Configuration ===================================================
fs = conf.fs;


%% ===== Computation =====================================================
a = cell(size(sos,1),1);
b = a;
for n=1:size(sos,1)
    if isoctave
        [bz,az] = bilinear(sos(n,1:3),sos(n,4:6),1/fs);
    else
        [bz,az] = bilinear(sos(n,1:3),sos(n,4:6),fs);
    end 
    b{n} = bz;
    a{n} = az;
end
