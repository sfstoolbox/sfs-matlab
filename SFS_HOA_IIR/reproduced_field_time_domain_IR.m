% plot the reproduced sound field for given driving function
% S.Spors, 24.8.2010

clear all;

% load loudspeaker driving functions
[d,fs]=wavread('IRs\iir_pw.wav');
d = [ zeros(512, size(d,2)); d; zeros(512, size(d,2)) ];

% parameters
c = 343;            % speed of sound

R=1.50;             % radius of loudspeaker array
nLS=56;             % number of loudspeakers

resolution = 400;   % spatial resolution


% generate spatial grid
X = linspace(-2, 2, resolution);
Y = linspace(-2, 2, resolution);
[x, y] = meshgrid(X, Y);

% get loudspeaker positions
[x0,n0] = LSpos_circ(0,0,R,nLS);
  

for t_0 = 1800;     %1012       % time(s) to evaluate (in samples)
    
    p = zeros(size(x));
    d_reshaped = zeros(size(x));
    
    for l = 1 : nLS
        
        % distance of secondary source to receiver position
        r = sqrt( (x-x0(1,l)).^2 + (y-x0(2,l)).^2 );
        % propagation duration from secondary source to receiver position
        t = (t_0/fs - r./c) .* fs + 1; % in samples

        % shift d appropriately and interpolate between samples to obtain
        % the amplitude at the receiver position
        d_reshaped = reshape( interp1( (1:size(d, 1)), d(:,l), t, 'linear'), ...
                                           resolution, resolution );
                                       
        % add to reproduced sound field (and assure causality)
        p = p + d_reshaped ./ r;
       
    end
    

    %plot field
    figure;
    imagesc(X, Y, p);
    draw_loudspeakers(x0,n0);
    
    axis square;
    turn_imagesc;
    colorbar;
    xlabel('x');
    ylabel('y');
    pause(0.1);

end


