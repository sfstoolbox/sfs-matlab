function varargout = sound_field_mono(X,Y,Z,x0,src,D,f,conf)
%SOUND_FIELD_MONO simulates a monofrequent sound field for the given driving
%signals and secondary sources
%
%   Usage: [P,x,y,z] = sound_field_mono(X,Y,Z,x0,src,D,f,[conf])
%
%   Input parameters:
%       X           - x-axis / m; single value or [xmin,xmax]
%       Y           - y-axis / m; single value or [ymin,ymax]
%       Z           - z-axis / m; single value or [zmin,zmax]
%       x0          - secondary sources [n x 6] / m
%       src         - source model for the secondary sources. This describes the
%                     Green's function, that is used for the modeling of the
%                     sound propagation. Valid models are:
%                       'ps' - point source
%                       'ls' - line source
%                       'pw' - plane wave
%       D           - driving signals for the secondary sources [m x n]
%       f           - monochromatic frequency / Hz
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       P           - Simulated sound field
%       x           - corresponding x axis / m
%       y           - corresponding y axis / m
%       z           - corresponding z axis / m
%
%   SOUND_FIELD_MONO(X,Y,Z,x0,src,D,f,conf) simulates a sound field
%   for the given secondary sources, driven by the corresponding driving
%   signals. The given source model src is applied by the corresponding Green's
%   function for the secondary sources. The simulation is done for one
%   frequency in the frequency domain, by calculating the integral for P with a
%   summation.
%   
%   To plot the result use plot_sound_field(P,x,y,z).
%
%   References:
%       G. Williams (1999) - "Fourier Acoustics", Academic Press
%
%   see also: plot_sound_field, sound_field_mono_wfs_25d

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
nargmin = 7;
nargmax = 8;
narginchk(nargmin,nargmax);
isargvector(X,Y,Z,D);
isargsecondarysource(x0);
isargpositivescalar(f);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end
if size(x0,1)~=length(D)
    error(['%s: The number of secondary sources (%i) and driving ', ...
        'signals (%i) does not correspond.'], ...
        upper(mfilename),size(x0,1),length(D));
end


%% ===== Configuration ==================================================
% Plotting result
useplot = conf.plot.useplot;
showprogress = conf.showprogress;


%% ===== Computation ====================================================
% Create a x-y-grid
[xx,yy,zz,x,y,z] = xyz_grid(X,Y,Z,conf);
% Check what are the active axes to create an empty sound field with the correct
% size
[~,x1,x2,x3]  = xyz_axes_selection(x,y,z);
% Initialize empty sound field
P = squeeze(zeros(length(x3),length(x2),length(x1)));
% Integration over secondary source positions
for ii = 1:size(x0,1)

    % progress bar
    if showprogress, progress_bar(ii,size(x0,1)); end

    % ====================================================================
    % Secondary source model G(x-x0,omega)
    % This is the model for the secondary sources we apply.
    % The exact function is given by the dimensionality of the problem, e.g. a
    % point source for 3D
    G = greens_function_mono(xx,yy,zz,x0(ii,1:3),src,f,conf);

    % ====================================================================
    % Integration
    %              /
    % P(x,omega) = | D(x0,omega) G(x-x0,omega) dx0
    %              /
    %
    % see: Williams1993 p. 36
    % x0(ii,7) is a weight for the single secondary sources which includes for
    % example a tapering window for WFS or a weighting of the sources for
    % integration on a sphere.
    P = P + D(ii) .* G .* x0(ii,7);

end

% === Scale signal (at xref) ===
P = norm_sound_field_at_xref(P,x,y,z,conf);

% return parameter
if nargout>0, varargout{1}=P; end
if nargout>1, varargout{2}=x; end
if nargout>2, varargout{3}=y; end
if nargout>3, varargout{4}=z; end

% ===== Plotting =========================================================
if nargout==0 || useplot
    plot_sound_field(P,x,y,z,x0,conf);
end
