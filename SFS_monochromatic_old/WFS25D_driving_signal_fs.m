% computes driving signals for reproduction of a point source
% [Verhejen]
%
% S.Spors, 26.7.2007

function [D] = WFS25D_driving_signal_fs(xs,k,x0,n0,xref,lssel,ns)


D=zeros(1,length(x0));

xd = bsxfun(@minus,x0',xs);
dr=sqrt( xd(:,1).^2 + xd(:,2).^2 );

cwin = dot(xd,[cos(n0); sin(n0)]',2);


% loudspeaker selection
if(nargin>6)
    swin = dot(xd,repmat([cos(ns) sin(ns)],length(n0),1),2);
else
    swin = dot(xd,repmat(-xs,length(n0),1),2);
end


% calculate driving function for active speakers
if(lssel)
    idx = find(swin<0);
else
    idx = 1:length(x0);
end


for n=idx
    D(n) = sqrt(k/(2*pi*j)) * cwin(n) .* 1./dr(n).^(3/2) .* exp(j*k*dr(n));
end