% NFC-HOA modal filter in SOS-DF II structure for a virtual point source
% S.Spors, 28.4.2011
function [hd] = HOA25D_modal_filter_pw(R,order,fs)


c=343;

% compute normalized roots/zeros of spherical Hankel function
B=zeros(1,order+2);
A=B;

for n=0:order
        B(n+1) = factorial(2*order-n)/(factorial(order-n)*factorial(n)*2^(order-n));
end
B=B(end:-1:1);
%A(2) = (-1)^order;
A(2) = 1;

% find zeros/roots
z=roots(B);
p=roots(A);


% compute SOS coefficients of modal driving function
[sos,g] = zp2sos(p,z*c/R,2,'down', 'none');

% transform coefficients
for n=1:size(sos,1)
    %[bz,az] = impinvar(sos(n,1:3),sos(n,4:6),fs);
    [bz,az] = bilinear(sos(n,1:3),sos(n,4:6),fs);
    if(length(bz)==2)
        sos(n,2:3)=bz;
        sos(n,5:6)=az;
    else
        sos(n,1:3)=bz;
        sos(n,4:6)=az;
    end
end


% realize FOS/SOS as DF-II structure
hd = dfilt.df2sos(sos); 
hd.ScaleValues(end)=2*(-1)^order;