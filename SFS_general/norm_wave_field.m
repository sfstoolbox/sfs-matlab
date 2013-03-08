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
isargmatrix(P);
isargvector(x,y,z);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
xref = position_vector(conf.xref);
xysamples = conf.xysamples;


%% ===== Computation =====================================================
% Get our active axis
[dimensions,x1,x2] = xyz_axes_selection(x,y,z);
if ~dimensions(1)
    xref(1) = xref(2);
    xref(2) = xref(3);
elseif ~dimensions(2)
    xref(2) = xref(3);
end
    
% Use the half of the x axis and xref
[a,x1idx] = find(x1>xref(1),1);
[a,x2idx] = find(x2>xref(2),1);
% abs(x(1)-x(end))/xysamples gives us the maximum distance between to samples.
% If abs(x(xidx)-xref(1)) is greater this indicates that we are out of our
% bounds
if isempty(x1idx) || abs(x1(x1idx)-xref(1))>2*abs(x1(1)-x1(end))/xysamples
    error('%s: your used conf.xref is out of your X boundaries',upper(mfilename));
end
if isempty(x2idx) || abs(x2(x2idx)-xref(2))>2*abs(x2(1)-x2(end))/xysamples
    error('%s: your used conf.xref is out of your Y boundaries',upper(mfilename));
end
% Scale signal to 1
P = 1*P/abs(P(x2idx,x1idx));
