% reproduction of point source using SDM
clear all;
%close all;

f     = 1000;
omega = 2*pi*f;
c     = 343;
yref = 2;

% position of point source
ys   =  1;
xs   =  0;

% Loudspeaker distance
dx = 0.01

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spatial_interval = [-20 20];
N = 4001; % delta_x = 0.01

x = linspace(spatial_interval(1), spatial_interval(2), N);
y = linspace(-1, 3, 400);

kxs = (2*pi) / dx;

% positive frequencies
kx = linspace(0, kxs/2, fix(N/2)+1);
kx(1) = kx(2);
% negative frequencies
kx = [fliplr( -kx(2:end-1) ), kx];

%[k_x_m y_m] = meshgrid(k_x, Y);
%y_m = abs(y_m);

% Indices for the ... and evanescent part
idxpr = ((abs(kx) <= omega/c));
idxev = ((abs(kx) > omega/c));

% ====== Spectrum of driving function ====================================
D_kx = zeros(1,length(kx));

% point source
D_kx(idxpr) = exp(1i*kx(idxpr)*xs) .* ...
    besselh(0,2,sqrt( (omega/c).^2 - kx(idxpr).^2 ) * (yref-ys) ) ./ ...
    besselh(0,2,sqrt( (omega/c).^2 - kx(idxpr).^2 ) * yref );
    
D_kx(idxev) = exp(1i*kx(idxev).*xs) .* ...
    besselk(0,sqrt( kx(idxev).^2 - (omega/c).^2 ) * (yref-ys) ) ./ ...
    besselk(0,sqrt( kx(idxev).^2 - (omega/c).^2 ) * yref );
    
% return

% %%%%%%%%%%%%%%%%%%%%%%%% bandlimitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B = 50; % min 2
% B_offset = 0;
% 
% %%%%% lowpass %%%%%
% central_frequency = round( length( D_kx )/2 ) - B_offset;
% 
% indices = [1 : central_frequency - B , ...
%             central_frequency + B : length( D_kx )];
% D_kx(indices) = 0;
% % D_kx = D_kx .* exp(i .* k_x .* 3);
% 
% % %%%%% highpass %%%%
% % indices = [1:B, N-(B-2):N]; 
% % G_tilde_kxy(:, indices) = 0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = ifftx( D_kx, [], 2);


%%%%%%%%%%%%%%%% discretization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  use only every 20th point in space (to simulate a loudspeaker spacing of
%  20cm
indexes = 1:length(D);
indexes = indexes( mod(indexes-1, 20)~=0 );
D( indexes )   =   0;

D_kx = fftx( D, [], 2 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% omnidirectional source
G_kx = zeros(length(y),length(kx));

[K,Y] = meshgrid(kx(idxpr),y);
G_kx(:,idxpr) = -1i/4 * ...
    besselh(0,2,sqrt( (omega/c).^2 - K.^2 ) .* abs(Y) );

[K,Y] = meshgrid(kx(idxev),y);
G_kx(:,idxev) = 1/(2*pi) * ...
    besselk(0,sqrt( K.^2 - (omega/c).^2 ) .* abs(Y) );


P_kx = repmat( D_kx, [size(G_kx, 1) 1] ) .* G_kx;
P    = ifftx( P_kx, [], 2);

% normalization
P    =    P ./ abs(P(end/2, end/2));

figure;
%GraphDefaults('paper');
imagesc(x, y, real(P), [-3 3])
colormap gray;
%imagesc(X, Y, db(abs(P)))
%imagesc(db(abs(G_kx)))
xlim([-2 2]);
xlabel('x');
ylabel('y');
%turn_imagesc; % Script by Jens?
axis square;



figure;
%GraphDefaults('paper');
plot(kx, 20*log10(abs(D_kx)) + 28, 'k') % + 28
xlim([-30 30]);
ylim([-50 10]);
xlabel('x');
axis square;
%GraphDefaults('paper');

% print -depsc2 ../../images/sound_fields_linear_continuous.eps
% print -depsc2 ../../images/sound_fields_linear_discrete.eps
% print -depsc2 ../../images/sound_fields_linear_discrete_1.eps
% print -depsc2 ../../images/sound_fields_linear_discrete_2.eps

% print -depsc2 ../../images/d_kx_point_source_1000Hz.eps
% print -depsc2 ../../images/d_kx_point_source_1000Hz_repetitions.eps
% print -depsc2 ../../images/d_kx_point_source_1000Hz_repetitions_bandlimited_1.eps
% print -depsc2 ../../images/d_kx_point_source_1000Hz_repetitions_bandlimited_2.eps
