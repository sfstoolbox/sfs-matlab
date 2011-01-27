% computes driving signals for reproduction of a plane wave
% using the 2.5-dimensional spectral division method
% CAUTION: works only for a linear/planar array
% S.Spors, 8.1.2010

function [D] = SDM25D_driving_signal_ps_approx(xs,k,x0,n0,xref,lssel)


D=zeros(1,length(x0));

% distance between source position and secondary sources
xd = bsxfun(@minus,x0',xs);
dr=sqrt( xd(:,1).^2 + xd(:,2).^2 );


for n=1:length(x0)
    D(n) = sqrt(xref(2)/(xref(2)-xs(2))) .* j*k .* xs(2) ./ dr(n) .* besselh(1,2,k*dr(n)); 
end