function [brs,delay] = ssr_brs_wfs(X,phi,xs,src,irs,conf)
%SSR_BRS_WFS generates a binaural room scanning (BRS) set for use with the
%SoundScape Renderer
%
%   Usage: [brs,delay] = ssr_brs_wfs(X,phi,xs,src,irs,conf)
%
%   Input parameters:
%       X       - listener position / m
%       phi     - azimuthal head orientation / rad
%       xs      - virtual source position [ys > Y0 => focused source] / m
%       src     - source type: 'pw' - plane wave
%                              'ps' - point source
%                              'fs' - focused source
%       irs     - impulse response data set for the secondary sources
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       brs     - conf.N x 2*nangles matrix containing all impulse responses (2
%                 channels) for every angles of the BRS set
%       delay   - delay added by driving function / s
%
%   SSR_BRS_WFS(X,phi,xs,src,irs,conf) prepares a BRS set for a virtual source
%   at xs for WFS and the given listener position. One way to use this BRS set
%   is using the SoundScapeRenderer (SSR), see http://spatialaudio.net/ssr/
%
%   See also: ir_generic, ir_wfs, driving_function_imp_wfs

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
isargposition(X);
isargxs(xs);
isargscalar(phi);
isargstruct(conf);


%% ===== Computation =====================================================
% Secondary sources
x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,xs,src);
x0 = secondary_source_tapering(x0,conf);
% Calculate driving function
[d,~,~,delay] = driving_function_imp_wfs(x0,xs,src,conf);
% Calculate brs set
brs = ssr_brs(X,phi,x0,d,irs,conf);
