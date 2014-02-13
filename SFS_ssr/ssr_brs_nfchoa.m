function brs = ssr_brs_nfchoa(X,phi,xs,src,irs,conf)
%SSR_BRS_NFCHOA generates a binaural room scanning (BRS) set for use with the
%SoundScape Renderer
%
%   Usage: brs = ssr_brs_nfchoa(X,phi,xs,src,irs,[conf])
%
%   Input parameters:
%       X       - listener position / m
%       phi     - listener direction [head orientation] / rad
%       xs      - virtual source position [ys > Y0 => focused source] / m
%       src     - source type: 'pw' - plane wave
%                              'ps' - point source
%                              'fs' - focused source
%       irs     - IR data set for the second sources
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       brs     - conf.N x 2*nangles matrix containing all impulse responses (2
%                 channels) for every angles of the BRS set
%
%   SSR_BRS_NFCHOA(X,phi,xs,src,irs,conf) prepares a BRS set for a virtual source
%   at position xs for a virtual loudspeaker array driven by nearfield
%   compensated higher order Ambisonics (NFC-HOA) and the given listener position.
%   One way to use this BRS set is using the SoundScapeRenderer (SSR), see
%   http://spatialaudio.net/ssr/
%
%   see also: ir_generic, ir_nfchoa, driving_function_imp_nfchoa

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
isargxs(xs);
isargscalar(phi);
check_irs(irs);
if nargin<nargmax
    conf = SFS_config;
end
isargstruct(conf);


%% ===== Computation =====================================================
% secondary sources
x0 = secondary_source_positions(conf);
% calculate driving function
d = driving_function_imp_nfchoa(x0,xs,src,conf);
% calculate brs set
brs = ssr_brs(X,phi,x0,d,irs,conf);
