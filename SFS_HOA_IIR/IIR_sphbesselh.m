% IIR structure for computation of a spherical Hankel function
% S.Spors, 3.9.2010

 
order=2;    %  order of spherical Hankel function
R=1.5;
c=343;
fs=44100;
f=linspace(0,20000,1000);
k=2*pi*f/c;
r=1.5;

% reference: spherical Hankel function
h1=sphbesselh(order,2,k*R);


% coefficient vectors for IIR filter 
B=zeros(1,order+2);
A=zeros(1,order+2);

for n=0:order
        B(n+1) = beta1(order,n);
        B(n+1) = B(n+1).* (R./c)^n;
end
B=B(end:-1:1);
A(1) = (R/c)^(n+1);

Z=A;
X=B;

TF1 = tf(B,A);
sysd1 = c2d(TF1,1/44100,'imp');

% calculate frequency response
h2=freqs(B,A,2*pi*f);

% multiply with remaining terms for comparison with reference
h2=-1i^n * h2 .* exp(-1i*2*pi*f*R/c);
%=========================================================================

%rect fct.
% t = 1:1:20000;
% d = 2:19998;
% Q = pulstran(t,d,'rectpuls');
% figure
% plot(t,Q)
% axis([-10 +350 -2 +2])
% 
% % function 1/(-ikr) for comparision
f=linspace(1,22050,1000);
k=2*pi.*f/c;
x=2./(-1i*k*r);
% 
% % numerator , denominator in s-domain
num = [0 0 0 0 0 0 0 2];
denum = (r/c) * [0 0 0 0 0 0 1 0];
% 
% % build transferfunction 1/s
TF = tf(num,denum);
[z,u] = freqs(num,denum,f);
sysd = c2d(TF,1/44100,'tustin');

z = z.*c/(2*pi*r);
% figure
% hold on
% plot(u,db(abs(z)),'b--'); % plot frequency response of z
% plot(f,db(abs(x)),'r:') % comparision with term 1/-ikr
% hold off
% 
% % transform from s to z domain via bilinear transf.
% % Integrator transf
[j,g] = bilinear(c*num,r*denum,fs);
[h,f1]=freqz(j,g,fs);
% figure
% hold on
% plot(f,db(abs(x)),'r:')
% plot(f1/pi*fs/2,db(abs(h)),'b--');
% hold off
% 
% % bessel fct transf.
[b,a] = bilinear(B,A,fs);
[h3,f3]=freqz(b,a,fs);

% 
% % impulse responses of the two filters
[o,p] = impz(j,g,[20000],fs);
figure
plot(p,o)

[u,i]= impz(a,b,[20000],fs);
figure
plot(i,u)

% convolution of impulse responses with cut of impulse resp. by rectangle window
I = conv(u,o);
% 
% %Fouriertransform of impulse repsonse
L = 300;                % Length of signal
NFFT = 2^nextpow2(L); % Next power of 2 from length of L
Y = fft(I,NFFT)/L;
f1 = fs/2*linspace(0,1,NFFT/2+1);

%Plot single-sided amplitude spectrum
figure
plot(f1,db(abs(Y(1:NFFT/2+1))))

%========================================================================




% show results
% figure
% subplot(2,1,1);
% 
% hold on
% plot(f,real(h1),'r-');
% plot(f,real(h2),'b--');
% hold off
% 
% subplot(2,1,2);
% 
% hold on
% plot(f,imag(h1),'r-');
% plot(f,imag(h2),'b--');
% hold off
% 
% 
% 
% discrete realization via bilinear transformation
% [b,a] = bilinear(B,A,fs);
% 
% [h3,f3]=freqz(b,a,fs);
% 
% [u,i]= impz(a,b,300,fs);
% figure
% plot(i,u)
% 
% 
% 
% figure
% hold on
% plot(f,db(abs(h1)),'r-');
% plot(f3/pi*fs/2,db(abs(h3)),'b--');
% hold off

%========================================================================


