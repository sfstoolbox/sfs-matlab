function varargout = wave_field_imp(X,Y,Z,x0,src,d,t,conf)
%WAVE_FIELD_IMP returns the wave field in time domain
%
%   Usage: [p,x,y,z] = wave_field_imp(X,Y,Z,x0,d,t,[conf])
%
%   Input options:
%       X           - x-axis / m; single value or [xmin,xmax]
%       Y           - y-axis / m; single value or [ymin,ymax]
%       Z           - z-axis / m; single value or [zmin,zmax]
%       x0          - positions of secondary sources / m
%       src         - source model of the Green's function. Valid models are:
%                       'ps' - point source
%                       'ls' - line source
%                       'pw' - plane wave
%       d           - driving function of secondary sources
%       t           - time / samples
%       conf        - optional configuration struct (see SFS_config)
%
%   Output options:
%       p           - simulated wave field
%       x           - corresponding x axis / m
%       y           - corresponding y axis / m
%       z           - corresponding z axis / m
%
%   WAVE_FIELD_IMP(X,Y,Z,x0,d,t,conf) computes the wave field synthesized by 
%   secondary sources driven by individual driving functions to the time t.
%   Point sources are applied as source models for the secondary sources.
%
%   To plot the result use:
%   conf.plot.usedb = 1;
%   plot_wavefield(p,x,y,z,conf);
%  
%   see also: wave_field_mono, plot_wavefield, greens_function_imp

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 7;
nargmax = 8;
narginchk(nargmin,nargmax);
isargvector(X,Y);
isargmatrix(x0,d);
isargchar(src);
isargscalar(t);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end
if size(x0,1)~=size(d,2)
    error(['%s: The number of secondary sources (%i) and driving ', ...
        'signals (%i) does not correspond.'], ...
        upper(mfilename),size(x0,1),size(d,2));
end

%% ===== Configuration ==================================================
% Plotting result
useplot = conf.plot.useplot;
% Speed of sound
c = conf.c;
% Sampling rate
fs = conf.fs;
% Debug mode
debug = conf.debug;
% Bandpass
usebandpass = conf.usebandpass;
bandpassflow = conf.bandpassflow;
bandpassfhigh = conf.bandpassfhigh;


%% ===== Computation =====================================================

% Spatial grid
[xx,yy,zz,x,y,z] = xyz_grid(X,Y,Z,conf);


% === Reshaping of the driving signal ===
%
% time reversal of driving function due to propagation of sound
% later parts of the driving function are emitted later by secondary
% sources
%
% ^      _     
% |     / \    driving function
% |    /   --
% | ---      --------
%  -----------------------------------> t
%
% ^            _
% |           / \  sound coming out of
% |         --   \ the secondary source
% | --------      ---
%  -----------------------------------> x
%   ^
%   position of secondary source
%
d = d(end:-1:1,:);
%
% Add additional zeros to the driving signal to ensure an amplitude of 0 in the
% whole listening area before and after the real driving signal.
% First get the maximum distance of the listening area and convert it into time
% samples
[dimensions,x1,x2,x3] = xyz_axes_selection(x,y,z); % get active axes
if all(dimensions)
    max_distance_in_samples = ...
        norm([x(1) y(1) z(1)]-[x(end) y(end) z(end)])/c * fs;
else
    max_distance_in_samples = norm([x1(1) x2(1)]-[x1(end) x2(end)])/c * fs;
end
% Append zeros at the beginning of the driving signal
d = [zeros(max_distance_in_samples,size(d,2)); d];
% correct time vector to work with inverted driving functions
% this will lead to a time point of t=0 for the starting of emitting the driving
% signal
t_inverted = t-size(d,1);
% Append zeros at the end of the driving signal
d = [d; zeros(max_distance_in_samples,size(d,2))];


% Initialize empty wave field (dependent on the axes we want)
p = squeeze(zeros(length(x3),length(x2),length(x1)));

% Integration over secondary sources
for ii = 1:size(x0,1)

    % Apply bandbass filter
    if usebandpass
        d(:,ii) = bandpass(d(:,ii),bandpassflow,bandpassfhigh,conf);
    end

    % ================================================================
    % Secondary source model: Greens function g3D(x,t)
    % distance of secondary source to receiver position
    [g,t_delta] = greens_function_imp(xx,yy,zz,x0(ii,1:3),src,t_inverted,conf);

    % Interpolate the driving function w.r.t. the propagation delay from
    % the secondary sources to a field point. The t returned from the Green's
    % function already inlcudes the desired time shift of the driving signal.
    % NOTE: the interpolation is required to account for the fractional
    % delay times from the loudspeakers to the field points
    ds = interp1(1:length(d(:,ii)),d(:,ii),t_delta,'spline');

    % ================================================================
    % Integration
    %          /
    % p(x,t) = | d(x0,t) g(x-x0,t) dx0
    %          /
    %
    % see: Spors2009, Williams1993 p. 36
    % x0(ii,7) is a weight for the single secondary sources which includes for
    % example a tapering window for WFS or a weighting of the sources for
    % integration on a sphere.
    p = p + ds .* g .* x0(ii,7);

end

% === Checking of wave field ===
check_wave_field(p,t);
% normalize field
p = norm_wave_field(p,conf);

% return parameter
if nargout>0 varargout{1}=p; end
if nargout>1 varargout{2}=x; end
if nargout>2 varargout{3}=y; end
if nargout>3 varargout{4}=z; end


% === Plotting ===
if nargout==0 || useplot
    plot_wavefield(p,x,y,z,x0,conf);
end

% some debug stuff
if debug
    figure; imagesc(db(abs(d))); title('driving functions'); caxis([-100 0]); colorbar;
end
