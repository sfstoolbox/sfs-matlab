function ir = ir_point_source(X,phi,xs,irs,conf)
%IR_POINT_SOURCE Generate a IR for a point source
%
%   Usage: ir = ir_point_source(X,phi,xs,irs,[conf])
%
%   Input parameters:
%       X       - listener position (m)
%       phi     - listener direction [head orientation] (rad)
%       xs      - source position (m)
%       irs     - IR data set for the second sources
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       ir      - Impulse response (nx2 matrix)
%
%   IR_POINT_SOURCE(X,phi,xs,irs,conf) calculates a impulse response for a
%   single loudspeaker at position xs and a listener located at X, looking
%   into direction phi. Whereby at phi = 0 the listener is looking in the
%   direction of the x-axis.
%
% see also: SFS_config, ir_wfs_25d, auralize_ir, brs_point_source

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input parameters ====================================
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
[X,xs] = position_vector(X,xs);
isargscalar(phi);
check_irs(irs);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Computation =====================================================
ir = ir_generic(X,phi,xs,1,irs,conf);
