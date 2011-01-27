% computes driving signals for 2.5D reproduction of a point source
% [Verhejen], [Spors, AES124]
%
% S.Spors, 26.7.2007

function [D] = WFS25D_driving_signal_ps(xs,k,x0,n0,xref,lssel)


D=zeros(1,length(x0));

xd = bsxfun(@minus,x0',xs);
dr=sqrt( xd(:,1).^2 + xd(:,2).^2 );

cwin = dot(xd,[cos(n0); sin(n0)]',2);

g0=1;

% calculate driving function for active speakers
if(lssel)
    idx = find(cwin>0);
else
    idx = 1:length(x0);
end


for n=idx    
    D(n) = g0 .* sqrt(k/(2*pi*j)) * cwin(n) .* 1./dr(n).^(3/2) .* exp(-j*k*dr(n));
end