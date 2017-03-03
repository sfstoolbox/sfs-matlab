function varargout = sound_field_mono(X,Y,Z,x0,src,D,f,conf)
%SOUND_FIELD_MONO simulates a monofrequent sound field for the given driving
%signals and secondary sources
%
%   Usage: [P,x,y,z] = sound_field_mono(X,Y,Z,x0,src,D,f,conf)
%
%   Input parameters:
%       X           - x-axis / m; single value or [xmin,xmax] or nD-array
%       Y           - y-axis / m; single value or [ymin,ymax] or nD-array
%       Z           - z-axis / m; single value or [zmin,zmax] or nD-array
%       x0          - secondary sources / m [nx7]
%       src         - source model for the secondary sources. This describes the
%                     Green's function, that is used for the modeling of the
%                     sound propagation. Valid models are:
%                       'ps' - point source
%                       'ls' - line source
%                       'pw' - plane wave
%       D           - driving signals for the secondary sources [mxn]
%       f           - monochromatic frequency / Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       P           - Simulated sound field
%       x           - corresponding x values / m
%       y           - corresponding y values / m
%       z           - corresponding z values / m
%
%   SOUND_FIELD_MONO(X,Y,Z,x0,src,D,f,conf) simulates a monochromatic sound
%   field for the given secondary sources, driven by the corresponding driving
%   signals. The given source model src is applied by the corresponding Green's
%   function for the secondary sources. The simulation is done for one
%   frequency in the frequency domain, by calculating the integral for P with a
%   summation.
%   For the input of X,Y,Z (DIM as a wildcard) :
%     * if DIM is given as single value, the respective dimension is
%     squeezed, so that dimensionality of the simulated sound field P is
%     decreased by one.
%     * if DIM is given as [dimmin, dimmax], a linear grid for the
%     respective dimension with a resolution defined in conf.resolution is
%     established
%     * if DIM is given as n-dimensional array, the other dimensions have
%     to be given as n-dimensional arrays of the same size or as a single value.
%     Each triple of X,Y,Z is interpreted as an evaluation point in an
%     customized grid.
%
%   To plot the result use:
%   plot_sound_field(P,X,Y,Z,conf);
%   or simple call the function without output argument:
%   sound_field_mono(X,Y,Z,x0,src,D,f,conf)
%
%   See also: plot_sound_field, sound_field_imp

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 8;
nargmax = 8;
narginchk(nargmin,nargmax);
isargnumeric(X,Y,Z);
isargvector(D);
isargsecondarysource(x0);
isargpositivescalar(f);
isargchar(src);
isargstruct(conf);
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
% Create a x-y-z-grid
[xx,yy,zz] = xyz_grid(X,Y,Z,conf);
[~,x1]  = xyz_axes_selection(xx,yy,zz); % get first non-singleton axis

% Initialize empty sound field
P = zeros(size(x1));

% Integration over secondary source positions
for ii = 1:size(x0,1)

    % progress bar
    if showprogress, progress_bar(ii,size(x0,1)); end

    % ====================================================================
    % Secondary source model G(x-x0,omega)
    % This is the model for the secondary sources we apply.
    % The exact function is given by the dimensionality of the problem, e.g. a
    % point source for 3D
    G = greens_function_mono(xx,yy,zz,x0(ii,1:6),src,f,conf);

    % ====================================================================
    % Integration
    %              /
    % P(x,omega) = | D(x0,omega) G(x-x0,omega) dx0
    %              /
    %
    % See http://sfstoolbox.org/#equation-single-layer
    %
    % x0(ii,7) is a weight for the single secondary sources which includes for
    % example a tapering window for WFS or a weighting of the sources for
    % integration on a sphere.
    P = P + D(ii) .* G .* x0(ii,7);

end

% return parameter
if nargout>0, varargout{1}=P; end
if nargout>1, varargout{2}=xx; end
if nargout>2, varargout{3}=yy; end
if nargout>3, varargout{4}=zz; end

% ===== Plotting =========================================================
if (nargout==0 || useplot)
    plot_sound_field(P,X,Y,Z,x0,conf);
end
