function brs = ssr_brs(X,phi,x0,d,irs,conf)
%SSR_BRS generates a binaural room scanning (BRS) set for use with the
%SoundScape Renderer
%
%   Usage: brs = ssr_brs(X,phi,x0,d,irs,conf)
%
%   Input parameters:
%       X       - listener position / m
%       phi     - azimuthal head orientation / rad
%       x0      - secondary sources
%       d       - corresponding driving signals
%       irs     - impulse response data set for the second sources
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       brs     - conf.N x 2*nangles matrix containing all impulse responses (2
%                 channels) for every angle of the BRS set
%
%   SSR_BRS(X,phi,x0,d,irs,conf) prepares a BRS set for the given secondary
%   sources and its driving signals for the given listener position.
%   One way to use this BRS set is using the SoundScapeRenderer (SSR), see
%   http://spatialaudio.net/ssr/
%
%   See also: ssr_brs_wfs, ssr_brs_nfchoa, ir_generic

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
isargsecondarysource(x0);
isargmatrix(d);
isargscalar(phi);
isargstruct(conf);


%% ===== Configuration ===================================================
N = conf.N;                       % Target length of BRIR impulse responses
angles = rad(conf.ir.brsangles);  % Angles for the BRIRs
showprogress = conf.showprogress; % Progress bar


%% ===== Computation =====================================================
nangles = length(angles);
% Initial values
brs = zeros(N,2*nangles);
% Generate a BRS set for all given angles
for ii = 1:nangles
    % Progress bar
    if showprogress, progress_bar(ii,nangles); end
    % Compute BRIR for the desired driving signals
    brs(:,(ii-1)*2+1:ii*2) = ...
        ir_generic(X,angles(ii)+phi,x0,d,irs,conf);
end
