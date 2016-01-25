function varargout = sound_field_imp_plane_wave(X,Y,Z,xs,t,conf)
%SOUND_FIELD_IMP_PLANE_WAVE simulates a sound field of a plane wave
%
%   Usage: [p,x,y,z] = sound_field_imp_plane_wave(X,Y,Z,xs,t,conf)
%
%   Input parameters:
%       X           - x-axis / m; single value or [xmin,xmax] or nD-array
%       Y           - y-axis / m; single value or [ymin,ymax] or nD-array
%       Z           - z-axis / m; single value or [zmin,zmax] or nD-array
%       xs          - direction of the plane wave / m
%       t           - time / samples
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       p           - Simulated sound field
%       x           - corresponding x values / m
%       y           - corresponding y values / m
%       z           - corresponding z values / m
%
%   SOUND_FIELD_IMP_PLANE_WAVE(X,Y,Z,xs,t,conf) simulates a sound
%   field of a plane wave going in the direction xs at time t.
%
%   To plot the result use:
%   plot_sound_field(p,X,Y,Z,conf);
%   or simple call the function without output argument:
%   sound_field_imp_plane_wave(X,Y,Z,xs,t,conf)
%   For plotting you may also consider to display the result in dB, by setting
%   the following configuration option before:
%   conf.plot.usedB = true;
%
%   See also: sound_field_imp, plot_sound_field, sound_field_mono_plane_wave

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
isargxs(xs);
isargstruct(conf);


%% ===== Computation ====================================================
% Disable the plotting of a source, because we have a plane wave
conf.plot.loudspeakers = 0;
[varargout{1:nargout}] = sound_field_imp(X,Y,Z,[xs 0 1 0 1],'pw',1,t,conf);
