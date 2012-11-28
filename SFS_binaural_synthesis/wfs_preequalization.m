function ir = wfs_preequalization(ir,conf)
%WFS_PREEQUALIZATION applies a pre-equalization filter for WFS
%
%   Usage: ir = wfs_prefilter(ir,[conf])
%
%   Input parameters:
%       ir      - IR to which the pre-equalization filter should be applied
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       ir      - IR with applied pre-equalization
%
%   WFS_PREEQUALIZATION(ir,conf) applies the pre-equalization filter for
%   Wave Field Synthesis to the given impulse response.
%
%   see also: wfs_prefilter, SFS_config, ir_wfs_25d

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
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
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
isargmatrix(ir);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
usehpre = conf.usehpre;


%% ===== Computation =====================================================
% Check if we should procide
if ~usehpre
    return;
end
% Get the filter
hpre = wfs_prefilter(conf);
% Apply the filter
for ii = 1:size(ir,2)
    ir(:,ii) = conv(hpre,ir(1:end-length(hpre)+1,ii));
end
