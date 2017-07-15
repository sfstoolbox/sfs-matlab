function brs = ssr_brs_point_source(X,phi,xs,irs,conf)
%SSR_BRS_POINT_SOURCE generates a binaural room scanning (BRS) set for use with
%the SoundScape Renderer
%
%   Usage: brs = ssr_brs_point_source(X,phi,xs,irs,conf)
%
%   Input parameters:
%       X       - listener position / m
%       phi     - azimuthal head orientation / rad
%       xs      - source position / m
%       irs     - impulse response data set for the second sources
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       brs     - conf.N x 2*nangles matrix containing all impulse responses (2
%                 channels) for every angle of the BRS set
%
%   SSR_BRS_POINT_SOURCE(X,phi,xs,irs,conf) prepares a BRS set for a reference
%   source (single point source) for the given listener position. One way to use
%   this BRS set is using the SoundScapeRenderer (SSR), see
%   http://spatialaudio.net/ssr/
%
%   See also: ir_generic, ir_point_source, get_ir

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
isargposition(X);
isargxs(xs);
isargscalar(phi);
isargstruct(conf);


%% ===== Computation =====================================================
brs = ssr_brs(X,phi,[xs 0 -1 0 1],1,irs,conf);
