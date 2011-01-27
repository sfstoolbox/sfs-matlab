% computes driving signals for reproduction of a line source
% S.Spors, 13.3.2008

function [D] = driving_signal_ls(xs,k,x0,n0,lssel)


D=zeros(1,length(x0));

xd = bsxfun(@minus,x0',xs);
dr=sqrt( xd(:,1).^2 + xd(:,2).^2 );

cwin = dot(xd,[cos(n0); sin(n0)]',2);


% calculate driving function for active speakers
if(lssel)
    idx = find(cwin>0);
else
    idx = 1:length(x0);
end


for n=idx
    D(n) = j*k * cwin(n) .* 1./dr(n) .*  besselh(1,2,k*dr(n));
end

