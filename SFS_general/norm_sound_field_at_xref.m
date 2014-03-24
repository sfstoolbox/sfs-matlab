function P = norm_sound_field_at_xref(P,x,y,z,conf)
%NORM_SOUND_FIELD_AT_XREF normalizes the sound field to 1 at xref
%
%   Usage: P = norm_sound_field_at_xref(P,x,y,z,[conf])
%
%   Input options:
%       P       - sound field
%       x,y,z   - vectors conatining the x-, y- and z-axis values / m
%       conf    - optional configuration struct (see SFS_config)
%
%   Output options:
%       P       - normalized sound field
%
%   NORM_SOUND_FIELD_AT_XREF(P,x,y,z,conf) normalizes the given sound field P to 1 at
%   the position conf.xref.
%
%   see also: norm_sound_field, sound_field_mono

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


%% ===== Checking of input parameters ====================================
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
isargnumeric(P);
isargvector(x,y,z);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
if ~conf.usenormalisation
    return;
end
xref = conf.xref;
resolution = conf.resolution;


%% ===== Computation =====================================================
% Get our active axis
[dimensions] = xyz_axes_selection(x,y,z);

% Use the half of the x axis and xref
if dimensions(1)
    xidx = find(x>xref(1),1);
    check_idx(xidx,x,xref(1),'X',resolution);
end
if dimensions(2)
    yidx = find(y>xref(2),1);
    check_idx(yidx,y,xref(2),'Y',resolution);
end
if dimensions(3)
    zidx = find(z>xref(3),1);
    check_idx(zidx,z,xref(3),'Z',resolution);
end

% Scale signal to 1
if all(dimensions)
    scale = abs(P(zidx,yidx,xidx));
elseif dimensions(1) && dimensions(2)
    scale = abs(P(yidx,xidx));
elseif dimensions(1) && dimensions(3)
    scale = abs(P(zidx,xidx));
elseif dimensions(2) && dimensions(3)
    scale = abs(P(zidx,yidx));
elseif dimensions(1)
    scale = abs(P(xidx));
elseif dimensions(2)
    scale = abs(P(yidx));
elseif dimensions(3)
    scale = abs(P(zidx));
else
    scale = 1;
end
P = P/scale;

end % of function


%% ===== Subfunctions ====================================================
function check_idx(idx,x,xref,str,resolution)
    % abs(x(1)-x(end))/resolution gives us the maximum distance between to samples.
    % If abs(x(xidx)-xref(1)) is greater this indicates that we are out of our
    % bounds
    if isempty(idx) || abs(x(idx)-xref)>2*abs(x(1)-x(end))/resolution
        error('%s: your used conf.xref is out of your %s boundaries', ...
            upper(mfilename),str);
    end
end
