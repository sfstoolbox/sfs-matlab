function Pm = circexp_mono_pw(npw,Nce,f,xq,conf)
%CIRCEXP_MONO_PW circular basis expansion of a mono-freqeunt plane wave
%
%   Usage: Pm = circexp_mono_pw(npw,Nce,xq,conf)
%
%   Input parameters:
%       npw     - propagation direction of plane wave [1 x 3]
%       Nce     - maximum order of circular basis expansion
%       f       - frequency of the monochromatic source / Hz
%       xq      - optional expansion center / m [1 x 3]
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       Pm      - regular circular expansion coefficients
%                 for m = 0:Nce, [1 x Nce+1]
%
%   See also: circexp_mono_ps

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


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 6;
narginchk(nargmin,nargmax);
isargcoord(xq,npw);
isargpositivescalar(Nce,f);
isargstruct(conf);


%% ===== Configuration ==================================================
c = conf.c;


%% ===== Computation =====================================================
[phipw, ~] = cart2pol(npw(1),npw(2));

m = (-Nce:Nce);
Pm = (-1i).^m.*exp(-1i.*m.*phipw);
Pm = Pm.*exp(-1i*xq*npw.'*2*pi*f./c);
