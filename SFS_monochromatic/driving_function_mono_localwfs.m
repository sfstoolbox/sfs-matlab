function [D, xv, x0] = driving_function_mono_localwfs(x0,xs,src,f,conf)
%DRIVING_FUNCTION_MONO_WFS returns the driving signal D for WFS
%
%   Usage: [D, xv, x0] = driving_function_mono_localwfs(x0,xs,src,f,conf)
%
%   Input parameters:
%       x0          - position and direction of the secondary source / m [nx6]
%       xs          - position of virtual source or direction of plane
%                     wave / m [1x3]
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs is the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       f           - frequency of the monochromatic source / Hz
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%       xv          - position and direction of the virtual secondary source / m [mx7]
%       x0          - position and direction of the secondary source / m [nx7]
%
%   see also: plot_sound_field, sound_field_mono_wfs

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
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
isargsecondarysource(x0);
isargxs(xs);
isargpositivescalar(f);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end

%% ===== Configuration ==================================================
virtualconf = conf;
virtualconf.secondary_sources.size = conf.virtual_secondary_sources.size;
virtualconf.secondary_sources.center = conf.virtual_secondary_sources.center;
virtualconf.secondary_sources.geometry = conf.virtual_secondary_sources.geometry;
virtualconf.secondary_sources.number = conf.virtual_secondary_sources.number;
%% ===== Computation ====================================================

% create virtual source array
xv = virtual_source_positions(x0,xs,src,conf);
% driving functions for virtual source array
Dv = driving_function_mono_wfs(xv,xs,src,f,virtualconf);

% select secondary sources for virtual secondary source array
selector = false(size(x0,1),1);
for idx=1:size(xv,1)
  [~, xdx] = secondary_source_selection(x0, xv(idx,1:6), 'fs');
  selector(xdx) = true;
end
x0(~selector,7) = 0;
% driving functions for real source array
D = driving_function_mono_wfs_vss(x0,xv,Dv,f,conf);

end
