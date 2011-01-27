% computes driving signals for reproduction of a plane wave
% with 2.5-dimensional higher-order Ambisonics
% S.Spors, 7.1.2010

function [D] = HOA25D_driving_signal_pw(al_pw,k,x0,R)

N=length(x0);
D=zeros(1,N);

al=atan2(x0(2,:),x0(1,:));


if mod(N,2)
   N21 =  floor(N/2);
   N22 = N21;
else
   N21 = N/2-1;
   N22 = N/2;
end


for n=0:N-1 % loop over all loudspeakers   
    for m = -N21 : N22
        D(n+1) = D(n+1) + 4.*pi .* ( j.^(-abs(m)) .* exp(-j.*m.*al_pw) ) ./ ( -j .* k .* sphbesselh(abs(m),2,k.*R) ) .*exp(j.*m.*al(n+1));
    end
end