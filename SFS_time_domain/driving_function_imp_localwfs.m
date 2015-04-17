function [d, x0, xv] = driving_function_imp_localwfs(x0,xs,src,conf)
%DRIVING_FUNCTION_IMP_LOCALWFS returns the driving signal d for local WFS
%
%   Usage: [D, xv, x0] = driving_function_mono_localwfs(x0,xs,src,conf)
%
%   Input parameters:
%       x0          - position and direction of the secondary source / m [nx6]
%       xs          - position of virtual source or direction of plane
%                     wave / m [1x3]
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs is the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'ls' - line source
%                         'fs' - focused source
%
%       f           - frequency of the monochromatic source / Hz
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%       x0          - position, direction, and weights of the real secondary
%                     sources / m [nx7]
%       xv          - position, direction, and weights of the virtual secondary
%                     sources / m [mx7]
%
%   References:
%       S. Spors (2010) - "Local Sound Field Synthesis by Virtual Secondary
%                          Sources", 40th AES
%
%   See also: plot_sound_field, sound_field_mono_wfs

%*****************************************************************************
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
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
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargsecondarysource(x0);
isargxs(xs);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
virtualconf = conf;
virtualconf.secondary_sources.size = conf.localsfs.vss.size;
virtualconf.secondary_sources.center = conf.localsfs.vss.center;
virtualconf.secondary_sources.geometry = conf.localsfs.vss.geometry;
virtualconf.secondary_sources.number = conf.localsfs.vss.number;
virtualconf.usetapwin = conf.localsfs.usetapwin;
virtualconf.tapwinlen = conf.localsfs.tapwinlen;
virtualconf.wfs = conf.localsfs.wfs;
method = conf.localsfs.method;


%% ===== Computation ====================================================
if strcmp('fs',src)
  error(['%s: %s is not a supported method source type! Try to use a point', ...
    ' source, if the source is inside the secondary source array but not', ...
    ' inside the virtual secondary source array'], upper(mfilename),src);
end

% Determine driving functions of virtual array with different sfs methods
switch method
  case 'wfs'
    % === Wave Field Synthesis ===
    % Create virtual source array
    xv = virtual_secondary_source_positions(x0,xs,src,conf);
    % Secondary_source_selection
    xv = secondary_source_selection(xv, xs, src);
    % Optional tapering
    xv = secondary_source_tapering(xv,virtualconf);
    % Optional amplitude correction
    % xv = secondary_source_amplitudecorrection(xv);
    % Driving functions for virtual source array
    dv = driving_function_imp_wfs(xv,xs,src,virtualconf);
  otherwise
    error('%s: %s is not a supported method for time domain localsfs!', ...
      upper(mfilename),method);
end

% Select secondary sources
x0 = secondary_source_selection(x0, xv(:,1:6), 'vss');
% Driving functions for real source array
d = driving_function_imp_wfs_vss(x0,xv,dv,conf);
end
