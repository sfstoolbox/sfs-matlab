% computes driving signals for reproduction of a plane wave 
% S.Spors, 26.7.2007

function [D] = HOA_driving_signal_pw(al_pw,k,x0,R)

N=length(x0);
D=zeros(1,N);

al=atan2(x0(2,:),x0(1,:));



for nu=-N/2:N/2
    D = D + j^(-nu) ./ besselh(nu,2,k.*R) .* exp(-j*nu*(al_pw-al));
end