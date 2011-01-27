% synesized wave field for 2.5D reproduction of a monochromatic 
% focused source using the spectral division method
% [Spors et al. ISCCSP 2010]

% S.Spors / 26.03.2010

f=1000;
c=343;
om=2*pi*f;

% position of focused source 
xs=0;
ys=1;

% loudspeaker distance
dx = 0.2;
%dx=0.01;

% init variables
kxrep=2*pi/dx;
Nrep=6;             % number of spectral repetitions

kxal=om/c;
%Nkx=1.5;            % factor by which kx is extended of kx=\omega/c criteria
Nkx = 4.5;

kx=linspace(-Nkx*kxal,Nkx*kxal,Nkx*2000);
y=linspace(0,3,300);
x=linspace(-2,2,300);

% Indexes for evanescent contributions and ...
%idxpr=find( abs(kx) <= (om/c) );
%idxev=find( abs(kx) > (om/c) );
idxpr=(( abs(kx) <= (om/c) ));
idxev=(( abs(kx) > (om/c) ));

% method to calculate driving function (only for non-aliased part)
withev=1;
yref = 2;



% ========== spectrum of driving function ========== 
Dpr=zeros(1,length(kx));
Dev=zeros(1,length(kx));


% non-aliased contributions
Dpr(idxpr) =  exp(-1j*kx(idxpr)*xs) .* ...
    besselh(0,2,sqrt((om/c)^2 - kx(idxpr).^2 )*(yref-ys)) ./ ...
    besselh(0,2,sqrt((om/c)^2 - kx(idxpr).^2 )*(yref));

if(withev)
    Dev(idxev) =  exp(-j*kx(idxev)*xs) .* ...
        besselk(0,sqrt(kx(idxev).^2 - (om/c).^2)*(yref-ys)) ./ ...
        besselk(0,sqrt(kx(idxev).^2 - (om/c).^2)*(yref));
end

%figure;plot(abs(Dpr'));
%figure;plot(abs(Dpr+Dev)');

% aliased contributions
Dpr_al=zeros(1,length(kx));
Dev_al=zeros(1,length(kx));

if(Nrep>0)

for n=-Nrep:Nrep
    
    if(n~=0)
    kxp=kx-n*kxrep;
    
    idx=(( abs(kxp) < (om/c) ));
    Dpr_al(idx) = Dpr_al(idx) + exp(-1j*kxp(idx)*xs) .* besselh(0,2,sqrt((om/c)^2 - kxp(idx).^2 )*(yref-ys)) ./ besselh(0,2,sqrt((om/c)^2 - kxp(idx).^2 )*(yref));
    
    if(withev)
        idx=(( abs(kxp) >= (om/c) ));
        Dev_al(idx) =  Dev_al(idx) + exp(-1j*kxp(idx)*xs) .* besselk(0,sqrt(kxp(idx).^2 - (om/c).^2)*(yref-ys)) ./ besselk(0,sqrt(kxp(idx).^2 - (om/c).^2)*(yref));
    end
    end
end

end


%  ==========  spectrum of secondary sources  ========== 
Gpr=zeros(length(kx),length(y));
Gev=zeros(length(kx),length(y));

[K,Y]=meshgrid(kx(idxpr),y);
Gpr(idxpr,:) = -1j/4.*besselh(0,2,sqrt( (om/c)^2 - K.^2 ).* Y)';


[K,Y]=meshgrid(kx(idxev),y);
Gev(idxev,:) = 1/(2*pi).*besselk(0,sqrt( K.^2 - (om/c)^2).* Y)';


%  ==========  reproduced wave field  ========== 
P=zeros(length(kx),length(y));
P_al=zeros(length(kx),length(y));

for n=1:length(y)
    % w/o aliasing
    P(idxpr,n) = Dpr(idxpr) .* Gpr(idxpr,n)';
    P(idxev,n) = Dev(idxev) .* Gev(idxev,n)';
    
    % aliasing
    P_al(idxpr,n) = Dpr_al(idxpr) .* Gpr(idxpr,n)';
    P_al(idxev,n) = Dev_al(idxev) .* Gev(idxev,n)';
end


%  ==========  inverse spatial Fourier transformation  ========== 
p=zeros(length(x),length(y));
p_al=zeros(length(x),length(y));

% Using Jens FFT
%p = ifftx(P', [], 2);
%p_al = ifftx(P_al', [], 2);
for n=1:length(x)
    for m=1:length(y)
        p(n,m) = sum ( P(:,m) .* exp(-1j*kx*x(n))' );
        p_al(n,m) = sum ( P_al(:,m) .* exp(-1j*kx*x(n))' );
    end
end

