f = 1000;

k = 2.*pi*f./343;
R = 1.5;
L = 56;
fs = 44100;

imp_resp = wavread('IRs\iir_pw.wav');

bin = round( f/44100 * size(imp_resp,1));


spec = ifft(imp_resp);

% for i = 1:56
%    spec(i,:) = spec(i,:).* exp(1i.*i.* ( +theta_pw - alpha_01(i) ) ); 
% end

x = linspace(-2, 2, 200);
y = linspace(-2, 2, 200);

P = zeros(200, 200);

for l = 1 : L
    
    D = spec(bin, l);
    alpha_0 = 2*pi/L * (l-1);
    P = P + D .*  point_source(x, y, R * cos(alpha_0), R * sin(alpha_0), k );
  
end
figure;
plot_wavefield;