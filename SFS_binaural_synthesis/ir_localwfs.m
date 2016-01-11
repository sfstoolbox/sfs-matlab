function [ir,x0] = ir_localwfs(X,phi,xs,src,irs,conf)
%IR_LOCALWFS generates a binaural simulation of local WFS
%
%   Usage: [ir,x0] = ir_localwfs(X,phi,xs,src,irs,conf)
%
%   Input parameters:
%       X       - listener position / m
%       phi     - listener direction [head orientation] / rad
%                 0 means the head is oriented towards the x-axis.
%       xs      - virtual source position / m
%       src     - source type: 'pw' -plane wave
%                              'ps' - point source
%       irs     - IR data set for the secondary sources
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       ir      - impulse response for the desired WFS array (nx2 matrix)
%       x0      - secondary sources / m
%
%   IR_LOCALWFS(X,phi,xs,src,irs,conf) calculates a binaural room impulse
%   response for a virtual source at xs for a virtual LOCAL WFS array and a
%   listener located at X.
%
%   See also: ssr_brs_wfs, ir_point_source, auralize_ir

%*****************************************************************************
% Copyright (c) 2010-2016 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2016 Institut fuer Nachrichtentechnik                   *
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
nargmin = 6;
nargmax = 6;
narginchk(nargmin,nargmax);
if conf.debug
    isargposition(X);
    isargxs(xs);
    isargscalar(phi);
    isargchar(src);
    isargstruct(conf);
end


%% ===== Computation =====================================================
% Get secondary sources
x0 = secondary_source_positions(conf);
% Get driving signals
[d, x0] = driving_function_imp_localwfs(x0,xs,src,conf);
% Generate the impulse response for WFS
ir = ir_generic(X,phi,x0,d,irs,conf);
