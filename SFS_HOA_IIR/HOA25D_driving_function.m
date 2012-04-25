% NFC-HOA driving function using IIR implementation
% S.Spors, 15.2.2011

clear all


% parameters
c=343;
fs=44100;

P=57;       % number of loudspeakers

r_ps=3;
theta_ps=0;

R=0.5;

N=1024;
N0=100;


% variables

if(isodd(P))
    order=(P-1)/2;    %  max order of spherical Hankel function
else
    order=floor((P+1)/2);
end
    
f=linspace(0,fs/2,N);
k=2*pi*f/c;
k(1)=k(2);


% compute impulse responses of modal filters
dm=zeros(order+1,N);
for n=1:order+1
    %df=HOA25D_modal_filter_ps(R,r_ps,n-1,fs);
    df=HOA25D_modal_filter_pw(R,n-1,fs);
    dm(n,:) = df.filter([zeros(1,N0) 1 zeros(1,N-1-N0)]);
end

% compute input signal for IFFT
d=zeros(2*order+1,N);

for n=-order:order
    d(n+order+1,:)=dm(abs(n)+1,:) .* exp(-1i*n*theta_ps);
end

% remove highest negative modal part for a even number of loudspeakers
if(iseven(P))
   d=d(2:end,:); 
end

if(1)
    % inverse Fourier transformation
    %d=ifftshift(d,1);
    d=circshift(d,[order+1 0]);
    d=(2*order+1)*ifft(d,[],1);
    d=d';
else
    % manual evaluation of Fourier series
    d_alpha_0 = 2*pi/(2*order+1);
    d2=d;
    
    for l = 0 : 2*order
        for m = -order : order    
            d(l+1,:) = d(l+1,:) + d2(m+order+1,:) .* exp(1i.*m.*l.*d_alpha_0);
        end
    end
    
    d=real(d)';
end


% plot results
t=1/fs*1000*(1:size(d,1));
d=d/max(abs(d(:)));

figure
%imagesc(1:(2*order+1),t,db(abs(d)));
imagesc(db(abs(d)));
turn_imagesc;
caxis([-100 0]);
tcolorbar('','dB');
xlabel('loudspeaker');
ylabel('time (ms)');
title('IIR filter representation');


figure
tmp=sum(d,2);
tmp=tmp/max(abs(tmp));
plot(t,tmp);
xlabel('time (ms)');
title('impulse response at center position');


if(0)
    %d2=driving_function_ambi(f, k, 2*order+1, order,theta_ps-pi,R);
    d2=driving_function_ambi_ps(f, k, 2*order+1, order,R,theta_ps,r_ps);
    d2=d2(1:2048,:);
    d2=d2/max(abs(d2(:)));
    
    t=1/fs*1000*(0:size(d2,1));

    figure
    imagesc(1:(2*order+1),t,db(abs(d2)));
    turn_imagesc;
    caxis([-100 0]);
    tcolorbar('','dB');
    xlabel('loudspeaker');
    ylabel('time (ms)');
    title('IIR filter representation');


    figure
    tmp=d2(2066-100:2066-100+1024-1,:);
    plot(90:120,d(90:120,:),'k-',90:120,tmp(90:120,:),'r-');
end


