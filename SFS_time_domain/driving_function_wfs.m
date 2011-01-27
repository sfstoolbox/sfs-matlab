function [d,ls_activity] = driving_function_wfs(X,Y,xs,ys,L,src,conf)
%DRIVING_FUNCTION_WFS

%% ===== Checking of input  parameters ==================================
nargmin = 6;
nargmax = 7;
error(nargchk(nargmin,nargmax,nargin));

if ~isnumeric(X) || ~isvector(X)
    error('%s: X has to be a vector!',upper(mfilename));
end
if ~isnumeric(Y) || ~isvector(Y)
    error('%s: Y has to be a vector!',upper(mfilename));
end
if ~isnumeric(xs) || ~isscalar(xs)
    error('%s: xs has to be a scalar!',upper(mfilename));
end
if ~isnumeric(ys) || ~isscalar(ys)
    error('%s: ys has to be a scalar!',upper(mfilename));
end
if ~isnumeric(L) || ~isscalar(L) || L<=0
    error('%s: L has to be a positive scalar!',upper(mfilename));
end
if ~ischar(src)
    error('%s: src has to be a string!',upper(mfilename));
end
if nargin<nargmax
    useconfig = true;
elseif ~isstruct(conf)
    error('%s: conf has to be a struct.',upper(mfilename));
else
    useconfig = false;
end


%% ===== Configuration ==================================================

% Load default configuration values
%if(useconfig)
if(useconfig)
    conf = SFS_config;
end
% Sampling rate
fs = conf.fs;
% Speed of sound
c = conf.c;
% Array position
X0 = conf.X0;
Y0 = conf.Y0;
% Array type
array = conf.array;
% Preequalization filter
usehpre = conf.usehpre;

%% ===== Computation =====================================================

maximum_delay = ceil(L/c * fs); % in samples

% WFS prefiltering
if (usehpre)
    prefilter = wfs_prefilter(conf);
else
    prefilter = 1;
end

% Get loudspeaker positions and directions
[x0,y0,phi] = secondary_source_positions(L,conf);
nLS = length(x0);
% Activity of secondary sources
ls_activity = secondary_source_selection(x0,y0,phi,xs,ys,src);

% Initialize driving function with maximum occuring delay plus some headroom
d = zeros(maximum_delay + 5*1024 + 2,nLS);

if strcmp('pw',src)
    % === Plane wave ===

    % Direction of plane wave
    nxs = xs / sqrt(xs^2+ys^2);
    nys = ys / sqrt(xs^2+ys^2);
    theta = -1*atan2(nxs,nys);
    % Direction of secondary sources
    nx0 = -sin(phi);
    ny0 = cos(phi);

    % delay in secs
    % NOTE: <n_pw,n(x0)> is the same as the cosinus between their angle
    %delay = L/c * (1-cos(phi(l) - theta));
    % FIXME: allign the frame positions for the different source types,
    % in order to plot them without the need to choose a different
    % frame.
    % FIXME: is it really needed to have different formulas for
    % different source types?
    if strcmp('circle',array)
        % Circular array
        % Delay for each single loudspeaker
        delay = L/2/c - L/2/c * (nxs*nx0 + nys*ny0);
    elseif strcmp('linear',array)
        % Linear array
        % TODO: express the sin term with nx0,ny0 in order to remove
        % theta.
        % Delay for each single loudspeaker to the wave front of the
        % plane wave.
        delay = L/c - abs(x0-x0(end))/c .* sin(phi-theta);
    elseif strcmp('box',array)
        to_be_implemented;
    else
        to_be_implemented;
    end
    % amplitude
    % FIXME: check this formula
    amplitude = cos(theta - phi);

elseif strcmp('ps',src)
    % === Point source ===
    % Distance between loudspeaker and virtual source
    r = sqrt((x0-xs).^2+(y0-ys).^2);
    % Delay and amplitude
    delay = r/c;
    amplitude = (ys-y0) .*  r.^(-2);
elseif strcmp('fs',src)
    % Focused source
    % Distance between loudspeaker and virtual source
    r = sqrt((x0-xs).^2+(y0-ys).^2);
    % Delay and amplitude (+max(r/c) ensures delay>=0)
    delay = L/c - r/c + max(r/c);
    amplitude = (ys-y0) .* r.^(-2);
else
    error('%s: %s is not a known source type.',upper(mfilename),src);
end

% delay in samples
delay = round(delay*fs);

% Iterate over loudspeaker and store driving signal
for l=1:nLS
    % Check if the desired secondary source is needed
    if ls_activity(l)>0
        % TODO: check if the multiplication with the prefilter is the same as a
        % convolution (dirac function).
        d(delay(l)+1:delay(l)+length(prefilter),l) = ...
            amplitude(l) .* prefilter;
    end
end

% put zeros around to have some headroom
d = [ zeros(512,nLS); d; zeros(512,nLS) ];

% normalize
d = d ./ max(abs(d(:)));
