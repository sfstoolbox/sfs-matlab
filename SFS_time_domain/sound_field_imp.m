function varargout = sound_field_imp(X,Y,Z,x0,src,d,t,conf)
%SOUND_FIELD_IMP returns the sound field in time domain
%
%   Usage: [p,x,y,z] = sound_field_imp(X,Y,Z,x0,src,d,t,conf)
%
%   Input options:
%       X           - x-axis / m; single value or [xmin,xmax] or nD-array
%       Y           - y-axis / m; single value or [ymin,ymax] or nD-array
%       Z           - z-axis / m; single value or [zmin,zmax] or nD-array
%       x0          - secondary sources / m [nx7]
%       src         - source model of the Green's function. Valid models are:
%                       'ps' - point source
%                       'ls' - line source
%                       'pw' - plane wave
%       d           - driving function of secondary sources
%       t           - time / s
%       conf        - configuration struct (see SFS_config)
%
%   Output options:
%       p           - simulated sound field
%       x           - corresponding x values / m
%       y           - corresponding y values / m
%       z           - corresponding z values / m
%
%   SOUND_FIELD_IMP(X,Y,Z,x0,src,d,t,conf) simulates a sound field for the
%   given secondary sources, driven by the corresponding driving signals. The
%   given source model src is applied by the corresponding Green's function
%   for the secondary sources. The simulation is done for the given time, by
%   calculating the integral for p with a summation.
%
%   To plot the result use:
%   plot_sound_field(p,X,Y,Z,conf);
%   or simple call the function without output argument:
%   sound_field_imp(X,Y,Z,x0,src,d,t,conf)
%   For plotting you may also consider to display the result in dB, by setting
%   the following configuration option before:
%   conf.plot.usedB = true;
%
%   See also: sound_field_mono, plot_sound_field, greens_function_imp

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
isargsecondarysource(x0);
isargchar(src);
isargmatrix(d);
isargscalar(t);
isargstruct(conf);
if size(x0,1)~=size(d,2)
    error(['%s: The number of secondary sources (%i) and driving ', ...
        'signals (%i) does not correspond.'], ...
        upper(mfilename),size(x0,1),size(d,2));
end


%% ===== Configuration ==================================================
c = conf.c;
fs = conf.fs;
L = conf.secondary_sources.size;
useplot = conf.plot.useplot;
debug = conf.debug;
showprogress = conf.showprogress;
usebandpass = conf.usebandpass;
bandpassflow = conf.bandpassflow;
bandpassfhigh = conf.bandpassfhigh;


%% ===== Computation =====================================================
% Spatial grid
[xx,yy,zz] = xyz_grid(X,Y,Z,conf);
[~,x1]  = xyz_axes_selection(xx,yy,zz); % get first non-singleton axis

% === Reshaping of the driving signal ===
%
% Time reversal of driving function due to propagation of sound later parts of
% the driving function are emitted later by secondary sources
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
% | x-------      ---
%  -----------------------------------> x
%   ^
%   position of secondary source
%
d = d(end:-1:1,:);

% Add additional zeros to the driving signal to ensure an amplitude of 0 in the
% whole listening area before and after the real driving signal.
%
% First get the maximum distance of the listening area and convert it into time
% samples, then compare it to the size 2*L of the secondary sources. If the size
% is larger use this for padding zeros.
max_listening_area = [max(xx(:)) max(yy(:)) max(zz(:))];
min_listening_area = [min(xx(:)) min(yy(:)) min(zz(:))];
max_distance = max(norm(max_listening_area-min_listening_area),2*L);
max_distance_in_samples = round(max_distance/c*fs);
% Append zeros at the beginning of the driving signal
d = [zeros(max_distance_in_samples,size(d,2)); d];
% Correct time vector to work with inverted driving functions, this will lead to
% a time point of t=0 for the starting of emitting the driving signal
t_inverted = t-size(d,1)/fs;
% Append zeros at the end of the driving signal
d = [d; zeros(max_distance_in_samples,size(d,2))];

% Initialize empty sound field (dependent on the axes we want)
p = zeros(size(x1));

% Apply bandpass filter
if usebandpass
    d = bandpass(d,bandpassflow,bandpassfhigh,conf);
end

% Integration over secondary sources
for ii = 1:size(x0,1)

    if showprogress, progress_bar(ii,size(x0,1)); end

    % ================================================================
    % Secondary source model: Greens function g3D(x,t)
    [g,t_delta] = greens_function_imp(xx,yy,zz,x0(ii,1:3),src,t_inverted,conf);

    % Interpolate the driving function w.r.t. the propagation delay from
    % the secondary sources to a field point. The t returned from the Green's
    % function already includes the desired time shift of the driving signal.
    % NOTE: the interpolation is required to account for the fractional
    % delay times from the loudspeakers to the field points
    ds = interp1(1:length(d(:,ii)),d(:,ii),t_delta*fs,'spline');

    % ================================================================
    % Integration
    %          /
    % p(x,t) = | d(x0,t) g(x-x0,t) dx0
    %          /
    %
    % See http://sfstoolbox.org/#equation-single-layer
    %
    % x0(ii,7) is a weight for the single secondary sources which includes for
    % example a tapering window for WFS or a weighting of the sources for
    % integration on a sphere.
    p = p + ds .* g .* x0(ii,7);

end

% === Checking of sound field ===
warning_if_zero(p,t);

% Return parameter
if nargout>0, varargout{1}=p; end
if nargout>1, varargout{2}=xx; end
if nargout>2, varargout{3}=yy; end
if nargout>3, varargout{4}=zz; end


%% ===== Plotting ========================================================
if nargout==0 || useplot
    plot_sound_field(p,X,Y,Z,x0,conf);
end

if debug
    figure; imagesc(db(abs(d))); title('driving functions'); caxis([-100 0]); colorbar;
end
