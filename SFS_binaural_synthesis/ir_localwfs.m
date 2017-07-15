function [ir,x0] = ir_localwfs(X,head_orientation,xs,src,sofa,conf)
%IR_LOCALWFS generates a binaural simulation of local WFS
%
%   Usage: [ir,x0] = ir_localwfs(X,head_orientation,xs,src,sofa,conf)
%
%   Input parameters:
%       X                - listener position / m
%       head_orientation - orientation of the listener with [phi theta] /
%                          (rad, rad)
%       xs               - virtual source position / m
%       src              - source type: 'pw' -plane wave
%                                       'ps' - point source
%       sofa             - impulse response data set for the secondary sources
%       conf             - configuration struct (see SFS_config)
%
%   Output parameters:
%       ir               - impulse responses for the desired WFS array
%                          (nx2 matrix)
%       x0               - secondary sources / m
%
%   IR_LOCALWFS(X,head_orientation,xs,src,sofa,conf) calculates a binaural room
%   impulse response for a virtual source at xs for a virtual LOCAL WFS array
%   and a listener located at X.
%
%   See also: ssr_brs_wfs, ir_point_source, auralize_ir

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
nargmin = 6;
nargmax = 6;
narginchk(nargmin,nargmax);
if conf.debug
    isargposition(X);
    isargxs(xs);
    isargvector(head_orientation);
    isargchar(src);
    isargstruct(conf);
end


%% ===== Computation =====================================================
% Get secondary sources
x0 = secondary_source_positions(conf);
% Get driving signals
[d, x0] = driving_function_imp_localwfs(x0,xs,src,conf);
% Generate the impulse response for WFS
ir = ir_generic(X,head_orientation,x0,d,sofa,conf);
