function D = driving_function_mono_localwfs_sbl(x0,xs,src,f,conf)
%DRIVING_FUNCTION_MONO_LOCALWFS_SBL returns the driving signal for local WFS
%using spatial bandwidth limitation
%
%   Usage: D = driving_function_mono_localwfs_sbl(x0,xs,src,f,conf)
%
%   Input parameters:
%       x0          - position and direction of the secondary source / m [nx7]
%       xs          - position of point source or direction of plane
%                     wave / m [1x3]
%       src         - source type of the virtual source
%                         'pw' - plane wave
%                         'ps' - point source
%       f           - frequency of the monochromatic source / Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving signals [mxn]
%
%   See also: sound_field_mono, sound_field_mono_localwfs_sbl

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
nargmin = 5;
nargmax = 5;
narginchk(nargmin,nargmax);
isargsecondarysource(x0);
isargxs(xs);
isargchar(src);
isargpositivescalar(f);
isargstruct(conf);


%% ===== Configuration ========================================================
xref = conf.xref;
N0 = size(x0,1);
% Resolution of plane wave decomposition
if isempty(conf.localwfs_sbl.Npw)
    Npw = 2*ceil(2*pi*0.9*f/conf.c*conf.secondary_sources.size/2);
else
    Npw = conf.localwfs_sbl.Npw;
end
% Maximum order of circular basis expansion of sound field
if isempty(conf.localwfs_sbl.order)
    Nce = nfchoa_order(N0,conf);
else
    Nce = conf.localwfs_sbl.order;
end


%% ===== Computation ==========================================================
% Circular expansion coefficients
switch src
case 'ps'
    Pm = circexp_mono_ps(xs,Nce,f,xref,conf);
case 'pw'
    Pm = circexp_mono_pw(xs,Nce,f,xref,conf);
otherwise
    error('%s: %s is not a known source type.',upper(mfilename),src);
end
% Modal window
wm = modal_weighting(Nce,conf);
Pm = bsxfun(@times,[wm(end:-1:2),wm],Pm);
% Plane wave decomposition
Ppwd = pwd_mono_circexp(Pm,Npw);
% Driving signal
D = driving_function_mono_wfs_pwd(x0,Ppwd,f,xref,conf);
