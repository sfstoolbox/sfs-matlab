% reproduction of point source using SDM
clear all;
%close all;

f     = 2000;
omega = 2*pi*f;
c     = 343;
y_ref = 2;

% position of point source
y_s   =  1;
x_s   =  0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spatial_interval = [-20 20];
N = 4001; % delta_x = 0.01

X = linspace(spatial_interval(1), spatial_interval(2), N);
Y = linspace(-1, 3, 400);

delta_x = ( spatial_interval(2) - spatial_interval(1) ) / N;
k_x_s = (2*pi) / delta_x;

% positive frequencies
k_x = linspace(0, k_x_s/2, N/2+1);
k_x(1) = k_x(2);

% negative frequencies
k_x = [fliplr( -k_x(2:end-1) ), k_x];

[k_x_m y_m] = meshgrid(k_x, Y);
y_m = abs(y_m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D_kx = zeros( size(k_x) );

% point source
D_kx( abs(k_x) <= omega/c ) = exp(i.*k_x( abs(k_x) <= omega/c ).*x_s) .* ...
    besselh( 0, 2, sqrt( (omega/c).^2 - k_x( abs(k_x) <= omega/c ).^2 ) * (y_ref-y_s) ) ./ ...
        besselh( 0, 2, sqrt( (omega/c).^2 - k_x( abs(k_x) <= omega/c ).^2 ) * y_ref );
    
D_kx( abs(k_x) > omega/c ) = exp(i.*k_x( abs(k_x)> omega/c ).*x_s) .* ...
    besselk(0, sqrt( k_x( abs(k_x) > omega/c ).^2 - (omega/c).^2 ) * (y_ref-y_s) ) ./ ...
       besselk(0, sqrt( k_x( abs(k_x) > omega/c ).^2 - (omega/c).^2 ) * y_ref );
    
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

figure; plot(D_kx)
figure; plot(real(D_kx))
figure; plot(imag(D_kx))

D = ifftx( D_kx, [], 2);


%%%%%%%%%%%%%%%% discretization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  use only every 20th point in space (to simulate a loudspeaker spacing of
%  20cm
%indexes = 1:length(D);
%indexes = indexes( mod(indexes-1, 20)~=0 );
%D( indexes )   =   0;

D_kx = fftx( D, [], 2 );

figure; plot(D_kx)
figure; plot(real(D_kx))
figure; plot(imag(D_kx))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% omnidirectional source
G_kx = zeros( size(y_m) );

G_kx( abs(k_x_m) <= omega/c ) = -i/4 * ...
    besselh(0, 2, sqrt( (omega/c).^2 - k_x_m( abs(k_x_m) <= omega/c ).^2 ) .* y_m( abs(k_x_m) <= omega/c ) );

G_kx( abs(k_x_m) > omega/c ) = 1/(2*pi) * ...
    besselk(0, sqrt( k_x_m( abs(k_x_m) > omega/c ).^2 - (omega/c).^2 ) .* y_m( abs(k_x_m) > omega/c ) );


P_kx = repmat( D_kx, [size(G_kx, 1) 1] ) .* G_kx;
P    = ifftx( P_kx, [], 2);

% normalization
P    =    P ./ abs(P(end/2, end/2));

figure;
%GraphDefaults('paper');
imagesc(X, Y, real(P), [-3 3])
colormap gray;
%imagesc(X, Y, db(abs(P)))
%imagesc(db(abs(G_kx)))
xlim([-2 2]);
xlabel('x');
ylabel('y');
%turn_imagesc;
axis square;



figure;
%GraphDefaults('paper');
plot(k_x, 20*log10(abs(D_kx)) + 28, 'k') % + 28
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
