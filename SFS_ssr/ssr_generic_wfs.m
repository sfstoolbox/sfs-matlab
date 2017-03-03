function ir = ssr_generic_wfs(xs,src,conf)
%SSR_GENRIC_WFS generates an impulse response for the generic renderer of the
%SoundScape Renderer
%
%   Usage: ir = ssr_generic_wfs(xs,src,conf)
%
%   Input parameters:
%       xs      - virtual source position / m
%       src     - source type: 'pw' -plane wave
%                              'ps' - point source
%                              'fs' - focused source
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       ir      - impulse response for the desired loudspeaker array
%
%   SSR_GENERIC_WFS(xs,src,conf) calculates an impulse response for a virtual
%   source at xs for the loudspeakers of a WFS array. Every loudspeaker of
%   the array is represented by one column in the impulse response.
%   For the generic renderer it is of importance to know what position the first
%   loudspeaker of your array will have. For example, the included circular
%   array has its first loudspeaker at phi=0deg which is on the x-axis. If you
%   have another setup you have to provide it with conf.secondary_sources.x0.
%
% See also: generic_nfchoa, brs_wfs, driving_function_imp_wfs

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
isargxs(xs);
isargchar(src);
isargstruct(conf);


%% ===== Configuration ==================================================
N = conf.N;


%% ===== Main ============================================================
% Secondary sources
x0 = secondary_source_positions(conf);
% Create empty impulse response for all secondary sources
ir = zeros(N,size(x0,1));
[x0,idx] = secondary_source_selection(x0,xs,src);
x0 = secondary_source_tapering(x0,conf);
% Driving signals for the active speakers
ir(:,idx) = driving_function_imp_wfs(x0,xs,src,conf);
