function ir = ssr_generic_wfs(xs,src,conf)
%SSR_GENRIC_WFS generates an impulse response for the generic renderer of the
%SoundScape Renderer
%
%   Usage: ir = ssr_generic_wfs(xs,src,[conf])
%
%   Input parameters:
%       xs      - virtual source position / m
%       src     - source type: 'pw' -plane wave
%                              'ps' - point source
%                              'fs' - focused source
%       conf    - optional configuration struct (see SFS_config)
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
% see also: generic_nfchoa, brs_wfs, driving_function_imp_wfs

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
nargmin = 2;
nargmax = 3;
narginchk(nargmin,nargmax);
isargxs(xs);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
end
isargstruct(conf);


%% ===== Configuration ==================================================
N = conf.N;


%% ===== Main ============================================================
% Secondary sources
x0 = secondary_source_positions(conf);
% create empty impulse response for all secondary sources
ir = zeros(N,size(x0,1));
[x0,idx] = secondary_source_selection(x0,xs,src);
x0 = secondary_source_tapering(x0,conf);
% driving signals for the active speakers
ir(:,idx) = driving_function_imp_wfs(x0,xs,src,conf);
