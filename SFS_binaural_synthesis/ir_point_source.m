function ir = ir_point_source(X,head_orientation,xs,sofa,conf)
%IR_POINT_SOURCE generates a binaural simulation of a point source
%
%   Usage: ir = ir_point_source(X,head_orientation,xs,sofa,conf)
%
%   Input parameters:
%       X                - listener position / m
%       head_orientation - orientation of the listener with [phi theta] /
%                          (rad, rad)
%       xs               - source position / m
%       sofa             - impulse response data set (sofa struct/file)
%       conf             - configuration struct (see SFS_config)
%
%   Output parameters:
%       ir               - impulse responses (nx2 matrix)
%
%   IR_POINT_SOURCE(X,head_orientation,xs,sofa,conf) calculates a impulse
%   response for a single loudspeaker at position xs and a listener located
%   at X, looking into head_orientation. Whereby at head_orientation = [0 0]
%   the listener is looking in the direction of the x-axis.
%
%   See also: ssr_brs_point_source, get_ir, ir_wfs, auralize_ir

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


%% ===== Checking of input parameters ====================================
nargmin = 5;
nargmax = 5;
narginchk(nargmin,nargmax);
isargposition(X);
isargxs(xs);
isargvector(head_orientation);
isargstruct(conf);


%% ===== Computation =====================================================
ir = ir_generic(X,head_orientation,[xs 0 -1 0 1],1,sofa,conf);
