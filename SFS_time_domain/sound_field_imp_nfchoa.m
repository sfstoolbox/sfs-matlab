function varargout = sound_field_imp_nfchoa(X,Y,Z,xs,src,t,conf)
%SOUND_FIELD_IMP_NFCHOA returns the sound field in time domain for NFC-HOA
%
%   Usage: [p,x,y,z,x0] = sound_field_imp_nfchoa(X,Y,Z,xs,src,t,[conf])
%
%   Input options:
%       X           - x-axis / m; single value or [xmin,xmax]
%       Y           - y-axis / m; single value or [ymin,ymax] 
%       Z           - z-axis / m; single value or [zmin,zmax]
%       xs          - position of point source / m
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs, ys are the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%       t           - time point t / samples
%       conf        - optional configuration struct (see SFS_config)
%
%   Output options:
%       p           - simulated sound field
%       x           - corresponding x axis / m
%       y           - corresponding y axis / m
%       z           - corresponding z axis / m
%       x0          - positions and directions of the secondary sources / m
%
%   SOUND_FIELD_IMP_NFCHOA(X,Y,Z,xs,src,t,conf) simulates a sound field of the
%   given source type (src) using a NFC-HOA driving function at the time point t.
%
%   To plot the result use:
%   conf.plot.usedb = 1;
%   plot_sound_field(p,x,y,z,x0,conf);
%
%   see also: driving_function_imp_nfchoa, sound_field_mono_nfchoa

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
nargmin = 6;
nargmax = 7;
narginchk(nargmin,nargmax);
isargvector(X,Y,Z);
isargxs(xs);
isargchar(src);
isargscalar(t);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
if strcmp('2D',conf.dimension)
    greens_function = 'ls';
else
    greens_function = 'ps';
end


%% ===== Computation =====================================================
% Get secondary sources
x0 = secondary_source_positions(conf);
% Calculate driving function
d = driving_function_imp_nfchoa(x0,xs,src,conf);
% Calculate sound field
[varargout{1:min(nargout,4)}] = ...
    sound_field_imp(X,Y,Z,x0,greens_function,d,t,conf);
% Return secondary sources if desired
if nargout==5, varargout{5}=x0; end
