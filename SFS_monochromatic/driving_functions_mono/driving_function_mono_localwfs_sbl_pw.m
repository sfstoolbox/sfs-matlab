function D = driving_function_mono_localwfs_sbl_pw(x0,nk,f,conf)
%DRIVING_FUNCTION_MONO_LOCALWFS_SBL_PW driving signal for a plane wave using
%local WFS with spatial bandwidth limitation
%
%   Usage: D = driving_function_mono_localwfs_sbl_pw(x0,nk,f,conf)
%
%   Input parameters:
%       x0          - position and direction of the secondary source / m [nx7]
%       nk          - propagation direction of plane wave / m [1x3]
%       f           - frequency of the monochromatic source / Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function [nx1]
%
%   See also: sound_field_mono_localwfs_sbl, driving_function_mono_localwfs_sbl
%
%   References:
%       Winter, Hahn, Spors (2017), "Time-Domain Realisation of Model-Based
%       Rendering for 2.5D Local Wave Field Synthesis Using Spatial
%       Bandwidth-Limitation", 25th European Signal Processing Conference
%       (EUSIPCO), pp. 718-722, https://bit.ly/2yMjlOw

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
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);


%% ===== Configuration ========================================================
N0 = size(x0,1);
xref = conf.xref; 
% Maximum order of circular basis expansion of sound field
if isempty(conf.localwfs_sbl.order)
    Nce = nfchoa_order(N0,conf);
else
    Nce = conf.localwfs_sbl.order;
end
% Resolution of plane wave decomposition
if isempty(conf.localwfs_sbl.Npw)
    Npw = 2*ceil(2*pi*0.9*f/conf.c*conf.secondary_sources.size/2);
else
    Npw = conf.localwfs_sbl.Npw;
end


%% ===== Computation ==========================================================
% Circular expansion coefficients, Winter et al. (2017), eq. (12)
Pm = circexp_mono_pw(nk,Nce,f,xref,conf);
% Modal window
wm = modal_weighting(Nce,conf);
Pm = bsxfun(@times,[wm(end:-1:2),wm],Pm);
% Plane wave decomposition, inverse FT of Winter et al. (2017), eq. (10)
Ppwd = pwd_mono_circexp(Pm,Npw); 
% Driving signal, Winter et al. (2017), eq. (8)
D = driving_function_mono_wfs_pwd(x0,Ppwd,f,xref,conf);
