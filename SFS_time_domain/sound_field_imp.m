function varargout = sound_field_imp(X,Y,Z,x0,src,d,t,conf)
%SOUND_FIELD_IMP returns the sound field in time domain
%
%   Usage: [p,x,y,z] = sound_field_imp(X,Y,Z,x0,src,d,t,conf)
%
%   Input options:
%       X           - x-axis / m; single value or [xmin,xmax] or nD-array
%       Y           - y-axis / m; single value or [ymin,ymax] or nD-array
%       Z           - z-axis / m; single value or [zmin,zmax] or nD-array
%       x0          - positions of secondary sources / m
%       src         - source model of the Green's function. Valid models are:
%                       'ps' - point source
%                       'ls' - line source
%                       'pw' - plane wave
%       d           - driving function of secondary sources
%       t           - time / samples
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
%   for the secondary sources. The simulation is done at one time sample, by
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
%   References:
%       H. Wierstorf (2014) - "Perceptual Assessment of Sound Field Synthesis",
%       PhD thesis, TU Berlin
%       G. Williams (1999) - "Fourier Acoustics", Academic Press
%
%   See also: sound_field_mono, plot_sound_field, greens_function_imp

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
nargmin = 8;
nargmax = 8;
narginchk(nargmin,nargmax);
isargnumeric(X,Y,Z);
isargmatrix(x0,d);
isargchar(src);
isargscalar(t);
isargstruct(conf);
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
% Progress bar
showprogress = conf.showprogress;
% Bandpass
usebandpass = conf.usebandpass;
bandpassflow = conf.bandpassflow;
bandpassfhigh = conf.bandpassfhigh;
L = conf.secondary_sources.size;


%% ===== Computation =====================================================
% Spatial grid
[xx,yy,zz] = xyz_grid(X,Y,Z,conf);
[~,x1]  = xyz_axes_selection(xx,yy,zz); % get first non-singleton axis

% === Reshaping of the driving signal ===
%
% Time reversal of driving function due to propagation of sound
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
% samples, than compare it to the size of the secondary sources. If the size is
% biger use this for padding zeros.
max_distance = norm( [max(xx(:)) max(yy(:)) max(zz(:))] - ...
    [min(xx(:)) min(yy(:)) min(zz(:))] );
max_distance_in_samples = round(max(max_distance/c*fs,2*L/c*fs));

% Append zeros at the beginning of the driving signal
d = [zeros(max_distance_in_samples,size(d,2)); d];
% Correct time vector to work with inverted driving functions
% this will lead to a time point of t=0 for the starting of emitting the driving
% signal
t_inverted = t-size(d,1);
% Append zeros at the end of the driving signal
d = [d; zeros(max_distance_in_samples,size(d,2))];

% Initialize empty sound field (dependent on the axes we want)
p = zeros(size(x1));

% Apply bandbass filter
if usebandpass
    d = bandpass(d,bandpassflow,bandpassfhigh,conf);
end

% Integration over secondary sources
for ii = 1:size(x0,1)

    % progress bar
    if showprogress, progress_bar(ii,size(x0,1)); end

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
    % see Wierstorf et al. (2015), eq.(#single:layer) or Williams (1993), p. 36
    % x0(ii,7) is a weight for the single secondary sources which includes for
    % example a tapering window for WFS or a weighting of the sources for
    % integration on a sphere.
    p = p + ds .* g .* x0(ii,7);

end

% === Checking of sound field ===
check_sound_field(p,t);

% Return parameter
if nargout>0, varargout{1}=p; end
if nargout>1, varargout{2}=xx; end
if nargout>2, varargout{3}=yy; end
if nargout>3, varargout{4}=zz; end


%% ===== Plotting ========================================================
if nargout==0 || useplot
    plot_sound_field(p,X,Y,Z,x0,conf);
end

% Some debug stuff
if debug
    figure; imagesc(db(abs(d))); title('driving functions'); caxis([-100 0]); colorbar;
end
