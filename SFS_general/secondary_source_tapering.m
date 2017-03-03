function x0 = secondary_source_tapering(x0,conf)
%SECONDARY_SOURCE_TAPERING applies a tapering window to the secondary sources
%
%   Usage: x0 = secondary_source_tapering(x0,conf)
%
%   Input options:
%       x0          - secondary sources / m
%       conf        - configuration struct (see SFS_config)
%
%   Output options:
%       x0          - secondary sources / m, containing the applied tapering
%                     window in its weights in x0(:,7)
%
%   SECONDARY_SOURCE_TAPERING(x0,conf) applies a tapering window (depending on
%   the conf.usetapwin and conf.tapwinlen settings) to the secondary sources.
%   It is applied to the weights stored in x0(:,7).
%
%   See also: secondary_source_positions, tapering_window

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
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);


%% ===== Calculation ====================================================
% Apply tapering window to secondary sources
x0(:,7) = x0(:,7) .* tapering_window(x0,conf);
