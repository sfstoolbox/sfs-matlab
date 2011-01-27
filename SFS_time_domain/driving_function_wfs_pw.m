function [ d ] = driving_function_wfs_pw( L, theta_pw, ...
                                                d_alpha_0, r0, c, fs )


%% ===== Configuration ==================================================

% Load default configuration values
%if(useconfig)
if(1)
    conf = SFS_config;
end
fs = conf.fs;
c = conf.c;


%% ===== Computation =====================================================

maximum_delay = ceil(r0/c * fs); % in samples

% WFS prefiltering
%prefilter = wavread('wfs_prefilter_100_1800_44100.wav');
prefilter = wfs_prefilter(conf);


% Get loudspeaker positions and directions
conf.array = 'circle';
L = 2*r0;
[LSpos,LSdir] = secondary_source_positions(L,conf);
x0 = LSpos(1,:);
y0 = LSpos(2,:);
zeta = LSdir;
nLS = length(x0)
% maximum occuring delay plus some headroom
d = zeros(maximum_delay + 5*1024 + 2,nLS);


% calculate driving functions
%for l = 1 : 56 % loop over loudspeakers
d_alpha_0 = 2*pi/nLS;

for l = 1:nLS

    alpha_0 = (l-1) * d_alpha_0;

    % secondary source selection
    if ( cos(zeta(l) - theta_pw + pi) < 0 )
        % continue with the next step in the loop
        continue;
    end

    % delay in secs
    delay =  r0/c * (1-cos(zeta(l) - theta_pw + pi));
    % delay in samples
    delay = round(delay*fs);

    % amplitude
    % FIXME: check this formula
    amplitude = cos(theta_pw - zeta(l)+ pi);

    d(delay+1:delay+length(prefilter), l) = amplitude .* prefilter;

end


% put zeros around to have some headroom
d = [ zeros(512,nLS); d; zeros(512,nLS) ];

% normalize
d = d ./ max(abs(d(:)));

end
