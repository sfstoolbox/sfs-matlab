function brs = ssr_brs(X,phi,x0,d,irs,conf)
%SSR_BRS generates a binaural room scanning (BRS) set for use with the
%SoundScape Renderer
%
%   Usage: brs = ssr_brs(X,phi,x0,d,irs,[conf])
%
%   Input parameters:
%       X       - listener position / m
%       phi     - listener direction [head orientation] / rad
%       x0      - secondary sources
%       d       - corresponding driving signals
%       irs     - IR data set for the second sources
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       brs     - conf.N x 2*nangles matrix containing all impulse responses (2
%                 channels) for every angles of the BRS set
%
%   SSR_BRS(X,phi,x0,d,irs,conf) prepares a BRS set for the given secondary
%   sources and its driving signals for the given listener position.
%   One way to use this BRS set is using the SoundScapeRenderer (SSR), see
%   http://spatialaudio.net/ssr/
%
%   see also: ssr_brs_wfs, ssr_brs_nfchoa, ir_generic

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 6;
narginchk(nargmin,nargmax);
isargposition(X);
isargsecondarysource(x0);
isargmatrix(d);
isargscalar(phi);
check_irs(irs);
if nargin<nargmax
    conf = SFS_config;
end
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
    % progress bar
    if showprogress, progress_bar(ii,nangles); end
    % compute BRIR for the desired driving signals
    brs(:,(ii-1)*2+1:ii*2) = ...
        ir_generic(X,angles(ii)+phi,x0,d,irs,conf);
end
