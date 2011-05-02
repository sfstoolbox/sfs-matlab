% IIR structure for computation of a spherical Hankel function
% S.Spors, 3.9.2010

 
order=4;    %  order of spherical Hankel function
R=1.5;
c=343;
fs=44100;
f=linspace(0,20000,1000);
k=2*pi*f/c;
r=1.5;

xs_rel = ( config.xps(1) - config.xref(1) );
ys_rel = ( config.xps(2) - config.xref(2) );
rs_rel = sqrt(  xs_rel.^2 + ys_rel.^2 ); %  get distance from point source

% reference: spherical Hankel function
h1=sphbesselh(order,2,k*rs_rel);

% coefficient vectors for IIR filter 
B=zeros(1,order+2);
A=zeros(1,order+2);

for n=0:order
        B(n+1) = beta1(order,n);
        B(n+1) = B(n+1).* (rs_rel./c)^n;
end
B=B(end:-1:1);
A(1) = (rs_rel/c)^(n+1);

% calculate frequency response
h2=freqs(B,A,2*pi*f);

% multiply with remaining terms for comparison with reference
h2=-1i^n * h2 .* exp(-1i*2*pi*f*rs_rel/c);
%=========================================================================

% show results
figure
subplot(2,1,1);
hold on
plot(f,real(h1),'r-');
plot(f,real(h2),'b--');
hold off

subplot(2,1,2);
hold on
plot(f,imag(h1),'r-');
plot(f,imag(h2),'b--');
hold off

% discrete realization via bilinear transformation
[b,a] = bilinear(B,A,fs);

[h3,f3]=freqz(b,a,fs);

% [u,i]= impz(a,b,300,fs);
% figure
% plot(i,u)

figure
hold on
plot(f,db(abs(h1)),'r-');
plot(f3/pi*fs/2,db(abs(h3)),'b--');
hold off

%========================================================================

% %Fouriertransform of impulse repsonse
% L = 300;                % Length of signal
% NFFT = 2^nextpow2(L); % Next power of 2 from length of L
% Y = fft(I,NFFT);
% f1 = fs/2*linspace(0,1,NFFT/2+1);
% 
% %Plot single-sided amplitude spectrum
% figure
% plot(f1,db(abs(Y(1:NFFT/2+1))))

%========================================================================
