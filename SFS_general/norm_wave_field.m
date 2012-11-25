function P = norm_wave_field(P,x,y,z,conf)
%NORM_WAVE_FIELD normalizes the wave field to 1 at xref
%
%   Usage: P = norm_wave_field(P,x,y,z,[conf])
%
%   Input options:
%       P       - wave field
%       x,y,z   - vectors conatining the x-, y- and z-axis values
%       conf    - optional configuration struct (see SFS_config)
%
%   Output options:
%       P       - normalized wave field
%
%   NORM_WAVE_FIELD(P,x,y,z,yref) normalizes the given wave field P to 1 at
%   the position conf.xref.
%
%   see also: wave_field_mono_wfs_25d

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
nargmin = 4;
nargmax = 5;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(P);
isargvector(x,y,z);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
xref = position_vector(conf.xref);


%% ===== Computation =====================================================
% Get our active axis
[x1,x2,xref] = active_dimensions(x,y,z,conf);
% Use the half of the x axis and xref
[a,x1idx] = find(x1>=xref(1),1);
[a,x2idx] = find(x2>=xref(2),1);
if isempty(x1idx) || abs(x1(x1idx)-xref(1))>0.1
    error('%s: your used conf.xref is out of your boundaries',upper(mfilename));
end
if isempty(x2idx) || abs(x2(x2idx)-xref(2))>0.1
    error('%s: your used conf.xref is out of your boundaries',upper(mfilename));
end
% Scale signal to 1
P = 1*P/abs(P(x2idx,x1idx));
