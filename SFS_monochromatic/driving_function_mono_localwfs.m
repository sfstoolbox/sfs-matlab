function [D, x0, xv, idx] = driving_function_mono_localwfs(x0,xs,src,f,conf)
%DRIVING_FUNCTION_MONO_LOCALWFS returns the driving signal D for local WFS
%
%   Usage: [D, xv, x0, idx] = driving_function_mono_localwfs(x0,xs,src,f,conf)
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
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%       x0          - position, direction, and weights of the real secondary
%                     sources / m [nx7]
%       xv          - position, direction, and weights of the virtual secondary
%                     sources / m [mx7]
%       idx         - index of the selected sources from the original x0
%                     matrix [mx1]
%
%   References:
%       S. Spors (2010) - "Local Sound Field Synthesis by Virtual Secondary
%                          Sources", 40th AES
%
%   See also: plot_sound_field, sound_field_mono_wfs

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
nargmin = 5;
nargmax = 5;
narginchk(nargmin,nargmax);
isargsecondarysource(x0);
isargxs(xs);
isargpositivescalar(f);
isargchar(src);
isargstruct(conf);


%% ===== Configuration ==================================================
virtualconf = conf;
virtualconf.usetapwin = conf.localsfs.usetapwin;
virtualconf.tapwinlen = conf.localsfs.tapwinlen;
virtualconf.secondary_sources.size = conf.localsfs.vss.size;
virtualconf.secondary_sources.center = conf.localsfs.vss.center;
virtualconf.secondary_sources.geometry = conf.localsfs.vss.geometry;
virtualconf.secondary_sources.number = conf.localsfs.vss.number;
method = conf.localsfs.method;


%% ===== Computation ====================================================

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
    % Driving functions for virtual source array
    Dv = driving_function_mono_wfs(xv,xs,src,f,virtualconf);
  case 'nfchoa'
    % === Near-Field-Compensated Higher Order Ambisonics ===
    % Create virtual source array
    xv = secondary_source_positions(virtualconf);
    % Driving functions for virtual source array
    Dv = driving_function_mono_nfchoa(xv,xs,src,f,virtualconf);
  otherwise
    error('%s: %s is not a supported method for localsfs!',upper(mfilename),method);
end

% Select secondary sources
[x0, idx] = secondary_source_selection(x0, xv(:,1:6), 'vss');
% Driving functions for real source array
D = driving_function_mono_wfs_vss(x0,xv,Dv,f,conf);
