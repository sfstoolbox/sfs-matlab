function [xx,yy,zz] = xyz_grid(X,Y,Z,conf)
%XYZ_GRID returns a xyz-grid for the listening area
%
%   Usage: [xx,yy,zz] = xyz_grid(X,Y,Z,conf)
%
%   Input parameters:
%       X        - x-axis / m; single value or [xmin,xmax]
%       Y        - y-axis / m; single value or [ymin,ymax]
%       Z        - z-axis / m; single value or [zmin,zmax]
%       conf     - configuration struct (see SFS_config)
%
%   Output parameters:
%       xx,yy,zz - matrices representing the xyz-grid / m
%
%   XYZ_GRID(X,Y,Z,conf) creates a xyz-grid to avoid a loop in the sound field
%   calculation for the whole listening area.
%
%   See also: xyz_axes_selection, is_dim_custom, sound_field_mono

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


%% ===== Checking input parameters =======================================
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargnumeric(X,Y,Z);
isargstruct(conf);


%% ===== Configuration ====================================================
resolution = conf.resolution;


%% ===== Computation =====================================================
dims = {X,Y,Z};

if any( is_dim_custom(X,Y,Z) )
  xx = X;
  yy = Y;
  zz = Z;
else
  % Check which dimensions will be non singleton
  dimensions = xyz_axes_selection(X,Y,Z);
  % Create xyz-axes
  xyz_axes = {X(1),Y(1),Z(1)};
  % create regular grid in each non-singleton dimension
  xyz_axes(dimensions) = cellfun( @(D) linspace(D(1),D(2),resolution).', ...
    dims(dimensions),'UniformOutput',false );
  % Create xyz-grid
  grids = xyz_axes;
  if sum(dimensions)>=2
    % create 2D/3D grid
    [grids{dimensions}] = meshgrid(xyz_axes{dimensions});
  end
  [xx,yy,zz] = grids{:};
end
