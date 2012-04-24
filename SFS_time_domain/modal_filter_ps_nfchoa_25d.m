% NFC-HOA modal filter in SOS-DF II structure for a virtual point source
% S.Spors, 28.4.2011
function [hd] = HOA25D_modal_filter_ps(R,r_ps,order,fs)


c=343;


% compute normalized roots/zeros of spherical Hankel function
B=zeros(1,order+2);
A=B;

for n=0:order
        B(n+1) = factorial(2*order-n)/(factorial(order-n)*factorial(n)*2^(order-n));
end
B=B(end:-1:1);
A(1) = 1;

% find zeros/roots
z=roots(B);
    

% compute SOS coefficients of modal driving function
if(1)
    [sos,g] = zp2sos(z*c/r_ps,z*c/R,1,'up', 'none'); 
else
    p=roots(A);
    [sos0,g] = zp2sos(z,p,1,'down', 'none');
    
    for n=1:size(sos0,1)
        sos(n,1) = sos0(n,1);
        sos(n,4) = sos0(n,1);
        
        sos(n,2) = c/r_ps*sos0(n,2);
        sos(n,3) = (c/r_ps)^2*sos0(n,3);
        
        sos(n,5) = c/R*sos0(n,2);
        sos(n,6) = (c/R)^2*sos0(n,3);
    end
end

% transform coefficients
for n=1:size(sos,1)
    %[bz,az] = impinvar(sos(n,1:3),sos(n,4:6),fs);
    [bz,az] = bilinear(sos(n,1:3),sos(n,4:6),fs,1000);
    %[bz,az] = ciim_sos(sos(n,1:3),sos(n,4:6),fs);
    
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
hd.ScaleValues(end)=1/(2*pi*R);