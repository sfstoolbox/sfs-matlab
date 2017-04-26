function [pm,delay_offset] = circexp_imp_pw(npw,Nce,xq,conf)
%CIRCEXP_IMP_PW calculates the circular basis expansion of a plane wave
%
%   Usage: [pm,delay_offset] = circexp_imp_pw(npw,Nce,xq,conf)
%
%   Input parameters:
%       npw     - propagation direction of plane wave [1 x 3]
%       Nce     - maximum order of circular basis expansion
%       xq      - optional expansion center / m [1 x 3]
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       pm            - regular circular expansion coefficients in time domain
%                       for m = 0:Nce, [conf.N x Nce+1]
%       delay_offset  - additional added delay, so you can correct it
%
%   CIRCEXP_IMP_PW(npw,Nce,xq,conf) returns the circular basis expansion of
%   a plane wave.
%
%   See also: circexp_imp_ps

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
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargcoord(xq,npw);
isargpositivescalar(Nce);
isargstruct(conf);


%% ===== Configuration ==================================================
N = conf.N;
c = conf.c;


%% ===== Computation =====================================================
[phipw, ~] = cart2pol(npw(1),npw(2));

pulse = dirac_imp();
% Compute impulse responses for each mode m
m = (0:Nce);
pm = (-1j).^m.*exp(-1j*phipw*m)*pulse;
pm = [pm; zeros(N-size(pm,1),Nce+1)];

delay_offset = -(xq*npw.')/c;
