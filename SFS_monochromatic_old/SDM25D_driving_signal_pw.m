% computes driving signals for reproduction of a plane wave
% using the 2.5-dimensional spectral division method
% CAUTION: works only for a linear/planar array
% S.Spors, 8.1.2010

function [D] = SDM25D_driving_signal_pw(al_pw,k,x0,n0,xref,lssel)


D=zeros(1,length(x0));

% compute driving function
kx=k*cos(al_pw);
ky=k*sin(al_pw);

for n=1:length(x0)
    D(n) = 4*j*exp(-j*ky*xref(2))./besselh(0,2,ky*xref(2)).*exp(-j*(kx*x0(1,n)+ky*x0(2,n)));
end