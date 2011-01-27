% computes driving signals for reproduction of a plane wave 
% S.Spors, 26.7.2007

function [D] = driving_signal_pw(al_pw,k,x0,n0,lssel)


D=zeros(1,length(x0));

% cos-window from directional gradient
cwin = cos(al_pw - n0);

% calculate driving function for active speakers
if(lssel)
    idx = find(cwin>=0);
else
    idx = 1:length(x0);
end

% secondary source eq
cwin = 2*i*k .* cwin;

% compute driving function
kx=k*cos(al_pw);
ky=k*sin(al_pw);

for n=idx
    D(n) = cwin(n) .* exp(-j*(kx*x0(1,n)+ky*x0(2,n)));
end

