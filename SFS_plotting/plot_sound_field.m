function plot_sound_field(P,X,Y,Z,x0,conf)
%PLOT_SOUND_FIELD plot the given sound field
%
%   Usage: plot_sound_field(P,X,Y,Z,[x0],conf)
%
%   Input parameters:
%       P           - matrix containing the sound field in the format P = P(y,x)
%       X           - x-axis / m; single value or [xmin,xmax] or nD-array
%       Y           - y-axis / m; single value or [ymin,ymax] or nD-array
%       Z           - z-axis / m; single value or [zmin,zmax] or nD-array
%       x,y,z       - vectors for the x-, y- and z-axis
%       x0          - matrix containing the secondary source positions to plot.
%                     Default: plot no secondary sources
%       conf        - configuration struct (see SFS_config)
%
%   PLOT_SOUND_FIELD(P,X,Y,Z,x0,conf) plots the sound field P in dependence of
%   the axes that are non-singleton. You have to provide the axes in the same
%   syntax as for the sound field functions. For a given set x0 of secondary
%   sources the secondary sources are added as dots or loudspeaker symbols
%   depending on your setting of conf.plot.realloudspeakers.
%
%   See also: sound_field_mono, sound_field_imp, draw_loudspeakers

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
nargmin = 5;
nargmax = 6;
narginchk(nargmin,nargmax);
isargnumeric(X,Y,Z,P);
if nargin<nargmax
    conf = x0;
    x0 = [];
end
if isempty(x0)
    conf.plot.loudspeakers = 0;
end
isargstruct(conf);


%% ===== Configuration ==================================================
% Tmp dir
tmpdir = conf.tmpdir;
% Plotting
p.usenormalisation = conf.plot.usenormalisation;
p.cmd = conf.plot.cmd;
p.usedb = conf.plot.usedb;
p.mode = conf.plot.mode;
p.size = conf.plot.size;
p.size_unit = conf.plot.size_unit;
p.caxis = conf.plot.caxis;
p.colormap = conf.plot.colormap;
p.loudspeakers = conf.plot.loudspeakers;
p.realloudspeakers = conf.plot.realloudspeakers;
p.lssize = conf.plot.lssize;
p.usefile = conf.plot.usefile;
p.file = conf.plot.file;


%% ===== Calculation =====================================================
% Handle the given axis and check which should be plotted
[dimensions,x1,x2,x3] = xyz_axes_selection(X,Y,Z);
% Check whether some of the axis have custom size
% NOTE: the grid will always be detected as custom, if x,y,z have more than two
% entries, which could also happen if you perform for example:
% >> [P,x,y,z,x0] = sound_field_mono_wfs(X,Y,Z,[0 -1 0],'pw',1000,conf);
% >> plot_sound_field(P,x,y,z,x0,conf)
% Instead you should do
% >> plot_sound_field(P,X,Y,Z,x0,conf)
% if X,Y,Z are a non-custom grid.
is_custom = is_dim_custom(X,Y,Z);

if any(is_custom) && sum(is_custom) < sum(dimensions)
  error(['%s: it is not possible to combine an axis with regular sampled ' ...
    'grid with an axis with custom grid']);
end

if ~any(dimensions)
    error(['%s: you have only one point in the sound field. ', ...
        'Omitting the plotting.'],upper(mfilename));
elseif ~dimensions(1) && ~dimensions(2)
    p.xlabel = 'z / m';
elseif ~dimensions(1)
    p.xlabel = 'y / m';
    p.ylabel = 'z / m';
elseif ~dimensions(2)
    p.xlabel = 'x / m';
    p.ylabel = 'z / m';
elseif ~dimensions(3)
    p.xlabel = 'x / m';
    p.ylabel = 'y / m';
else
    p.xlabel = 'x / m';
    p.ylabel = 'y / m';
    p.zlabel = 'z / m';
end

% Remove Inf values from sound field
P(isinf(P)) = NaN;

% Normalisation
if p.usenormalisation
    P = norm_sound_field(P,conf);
end

if p.usedb
    % If we have only zeros in the sound field set the field to eps to avoid
    % problems with log(0).
    if max(abs(P(:)))==0
        P(:) = eps;
    end
    % Change default colormap to yellowred
    if strcmp('default',p.colormap)
        conf.plot.colormap = 'yellowred';
    end
    % Calculate sound pressure level in dB
    P = 20*log10(abs(P));
    if p.usenormalisation
        P = P - max(P(:)); % ensure 0 dB max for dB plots
    end
else
    P = real(P);
end

% Check if we should plot loudspeakers symbols
if p.realloudspeakers && size(x0,1)>1000
    warning(['%s: the given number of loudspeaker is >1000. ',...
            'Switching back to non loudspeaker symbols.'],upper(mfilename));
    p.realloudspeakers = 0;
end

% Set the color bar axis to default values if not given otherwise
if p.caxis, else
    if p.usedb
        p.caxis = [-45,0];
    else
        p.caxis = [-1,1];
    end
end

%% ===== Plotting ========================================================

% Create a new figure
figure;
% Set size
figsize(p.size(1),p.size(2),p.size_unit);

switch sum(dimensions)
  case 1
    % === 1D Plot ====
    if any(is_custom)  % custom grid
        plot(x1,P,'o','MarkerFaceColor','blue','MarkerSize',3)
    else % regular grid
        x1 = linspace(x1(1),x1(2),conf.resolution);
        plot(x1,P);
    end
    xlabel(p.xlabel);
    if p.usedb
        ylabel('Amplitude / dB');
    else
        ylabel('Amplitude');
    end
  case 2
    % === 2D Plot ====
    if any(is_custom)  % custom grid
        point_size = [];
        if isoctave, point_size = 2; end % fix scatter plot point size for Octave
        scatter(x1(:),x2(:),point_size,limit_colors(P,p.caxis),'filled');
    else  % regular grid
        imagesc(x1,x2,P,p.caxis);
    end
    % Add color bar
    set_colorbar(conf);
    % Set the y direction in normal mode (imagesc uses the reverse mode by
    % default)
    turn_imagesc;
    % Set the axis to use the same amount of space for the same length (m)
    axis image;
    % Labels etc. for the plot
    xlabel(p.xlabel);
    ylabel(p.ylabel);
    % Add loudspeaker to the plot
    if p.loudspeakers % && dimensions(1) && dimensions(2)
        hold on;
        draw_loudspeakers(x0,dimensions,conf);
        hold off;
    end
  case 3
    if ~any(is_custom)
      [x1,x2,x3] = xyz_grid(X,Y,Z,conf);
    end
    % === 3D Plot ====
    point_size = [];
    if isoctave, point_size = 5; end % fix scatter plot point size for Octave
    scatter3(x1(:),x2(:),x3(:),point_size,limit_colors(P,p.caxis),'filled');
    % Fix perspective of plot in Matlab (http://bit.ly/1LISaz9)
    set(gcf,'renderer','opengl');
    % Add color bar
    set_colorbar(conf);
    % Set the axis to use the same amount of space for the same length (m)
    axis image;
    % Labels etc. for the plot
    xlabel(p.xlabel);
    ylabel(p.ylabel);
    zlabel(p.zlabel);
end

% Save as file
if ~isempty(p.file) && strcmp('png',p.file(end-2:end))
    print_png(p.file);
elseif ~isempty(p.file) && strcmp('eps',p.file(end-2:end))
    print_eps(p.file);
end
end


%% ===== Helper functions ================================================
function P = limit_colors(P,caxis)
    % Limits the number of different values in P to 64.
    % This speeds up plotting with the scatter function in octave, see:
    % http://savannah.gnu.org/bugs/?40663
    %
    number_of_colors = 64;
    % First apply the caxis values
    P = min(caxis(2),max(caxis(1),P(:)));
    % Transform to [0...64] and transform to integer
    P = fix(number_of_colors/abs(max(P(:))-min(P(:))) .* (P+abs(min(P(:)))));
    % Transform back to [caxis(1)...caxis(2)]
    P = abs(caxis(2)-caxis(1))/number_of_colors .* P + caxis(1);
end
