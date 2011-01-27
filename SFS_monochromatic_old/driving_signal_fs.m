% computes driving signals for reproduction of a focused point source
% S.Spors, 26.7.2007

function [D] = driving_signal_fs(xs,k,x0,n0,lssel,nfs)


D=zeros(1,length(x0));


dr=sqrt( (x0(1,:)-xs(1)).^2 + (x0(2,:)-xs(2)).^2 );
al=atan2(-x0(2,:)+xs(2),-x0(1,:)+xs(1));

% cos-window from directional gradient
cwin = cos(al - n0);

% calculate driving function for active speakers
if(lssel)
    cnfs=[cos(nfs);sin(nfs)];
    for n=1:length(x0)
        dp(:,n)=cnfs'*(xs'-x0(:,n));
    end
    idx = find(dp > 0);
else
    idx = 1:length(x0);
end


for n=idx
    D(n) = cwin(n) .* 1./dr(n) .* exp(j*k*dr(n));
%   D(n) = dr(n) .* exp(j*k*dr(n));
end