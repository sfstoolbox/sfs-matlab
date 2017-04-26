function [d,delay_offset] = driving_function_imp_localwfs_sbl_pw(x0,nk,conf)

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


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 3;
narginchk(nargmin,nargmax);


%% ===== Configuration ========================================================
N0 = size(x0,1);
xref = conf.xref; 
fs = conf.fs;
% maximum order of circular basis expansion of sound field
if isempty(conf.localsfs.sbl.order)
    Nce = nfchoa_order(N0,conf);
else
    Nce = conf.localsfs.sbl.order;
end
% resolution of plane wave decomposition
if isempty(conf.localsfs.sbl.Npw)
    Npw = 2*ceil(2*pi*0.9*fs/conf.c*conf.secondary_sources.size/2);
else
    Npw = conf.localsfs.sbl.Npw;
end

wfsconf = conf;
wfsconf.wfs = conf.localsfs.wfs;


%% ===== Computation ==========================================================
% circular expansion coefficients
[pm,delay_circexp] = circexp_imp_pw(nk,Nce,xref,conf);
% plane wave decomposition
ppwd = pwd_imp_circexp(pm,Npw);
% driving signal
[d,delay_lwfs] = driving_function_imp_wfs_pwd(x0,ppwd,xref,wfsconf);
% delay
delay_offset = delay_lwfs + delay_circexp;
