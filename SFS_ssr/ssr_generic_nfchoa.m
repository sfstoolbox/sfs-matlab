function ir = ssr_generic_nfchoa(xs,src,conf)
%SSR_GENERIC_NFCHOA generate an impulse response for the generic renderer of the
%SoundScape Renderer
%
%   Usage: ir = ssr_generic_nfchoa(xs,src,conf)
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
%   GENERIC_NFCHOA(xs,src,conf) calculates an impulse response for a virtual
%   source at xs for the loudspeakers of a NFC-HOA array. Every loudspeaker of
%   the array is represented by one column in the impulse response.
%   For the generic renderer it is of importance to know what position the first
%   loudspeaker of your array will have. For example, the included circular
%   array has its first loudspeaker at phi=0deg which is on the x-axis. If you
%   have another setup you have to provide it with conf.secondary_sources.x0.
%
%
%   See also: generic_wfs, brs_nfchoa, driving_function_imp_nfchoa

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


%% ===== Main ============================================================
% Secondary sources
x0 = secondary_source_positions(conf);
% Driving signals for the secondary sources
ir = driving_function_imp_nfchoa(x0,xs,src,conf);
