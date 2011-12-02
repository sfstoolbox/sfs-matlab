function [b,a]=ciim_sos(bi,ai,fs)

T=1/fs;

% second order section, complex conjugate zeros/poles
if(bi(3)~=0 && ai(3)~=0)
    z = roots(bi);
    p = roots(ai);

    alb = real(z(1));
    omb = imag(z(1));

    ala = real(p(1));
    oma = real(p(1));


    b(1) = ala - alb + 1/T;
    b(2) = ( ((ala-alb)^2-oma^2+omb^2)/oma * sin(oma*T) - 2/T * cos(oma*T))*exp(ala*T);
    b(3) = -(ala-alb-1/T)*exp(2*ala*T);

    a(1) = 1;
    a(2) = -2*cos(oma*T)*exp(ala*T);
    a(3) = exp(2*ala*T);

    b=T*b;
end

% first order section real zeros/poles
if(bi(3)==0 && ai(3)==0)

%     alb = -bi(2);
%     ala = -ai(2);
%     
%     b(1) = T/2 * (ala-alb)+1;
%     b(2) = (T/2*(ala-alb)-1)*exp(ala*T);
%     
%     a(1) = 1;
%     a(2) = - exp(ala*T);

    [b,a] = bilinear(bi,ai,fs);

end