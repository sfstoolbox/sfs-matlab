function [ d ] = driving_function_wfs_sw( L, x_s, y_s, ...
                                                d_alpha_0, r_0, c, fs )


r_s     = sqrt( x_s^2 + y_s^2 );
%alpha_s = atan2(y_s, x_s);


maximum_delay = ceil( (r_s+r_0) / c * fs); % in samples

% maximum occuring delay plus some headroom
d = zeros(maximum_delay + 5 * 1024 + 2 , L);

% WFS prefiltering
prefilter = wavread('wfs_prefilter_100_1800_44100.wav');

% calculate driving functions
for l = 1 : 56 % loop over loudspeakers

    alpha_0 = (l-1) * d_alpha_0;

    x_0 = r_0 * cos(alpha_0);
    y_0 = r_0 * sin(alpha_0);

    % vector between source position and ls position
    delta_r     = sqrt((x_0-x_s).^2 + (y_0-y_s).^2);
    delta_alpha = atan2(y_0-y_s, x_0-x_s);

    % secondary source selection
    if ( cos(alpha_0 - delta_alpha + pi) < 0 )
        continue;
    end

    % delay in secs
    delay =  delta_r / c;
    % delay in samples
    delay = round( delay * fs );

    amplitude = 1/delta_r * cos(alpha_0 - delta_alpha + pi);

    d(delay + 1 : delay + length(prefilter), l) = amplitude .* prefilter;

end


% put zeros around to have some headroom
d = [ zeros(512, L); d; zeros(512, L) ];

% normalize
d = d ./ max(abs(d(:)));

end
