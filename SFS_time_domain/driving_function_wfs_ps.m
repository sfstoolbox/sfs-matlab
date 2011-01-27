function [ d ] = driving_function_wfs_pw( L, xs, ys, d_alpha_0, r_0, c, fs )


maximum_delay = ceil(r_0 / c * fs); % in samples

% maximum occuring delay plus some headroom 
d = zeros(maximum_delay + 5 * 1024 + 2 , L);

% WFS prefiltering
prefilter = wavread('wfs_prefilter_100_1800_44100.wav');

% calculate driving functions
for l = 1 : 56 % loop over loudspeakers
    
    alpha_0 = (l-1) * d_alpha_0;
    
    % secondary source selection
    if ( cos(alpha_0 - theta_pw + pi) < 0 ) 
        continue;
    end

    % delay in secs
    delay =  r_0 / c * (1-cos(alpha_0 - theta_pw + pi));
    % delay in samples
    delay = round( delay * fs );
    
    amplitude = cos(theta_pw - alpha_0 + pi);
         
    d(delay + 1 : delay + length(prefilter), l) = amplitude .* prefilter;
       
end
  

% put zeros around to have some headroom
d = [ zeros(512, L); d; zeros(512, L) ];

% normalize
d = d ./ max(abs(d(:)));

end
