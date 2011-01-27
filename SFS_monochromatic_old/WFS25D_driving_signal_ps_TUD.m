% computes driving signals for 2.5D reproduction of a point source
% [Verhejen, Phd thesis, TUD], [Spors, AES124, Eq.(33)], [Opperschall, DA]
% CAUTION: only valid for linear array parallel to x-axis
%
% S.Spors, 4.2.2010

function [D] = WFS25D_driving_signal_ps_TUD(xs,k,x0,n0,xref,lssel)


D=zeros(1,length(x0));

% distance between source position and secondary sources
xd = bsxfun(@minus,x0',xs);
dr=sqrt( xd(:,1).^2 + xd(:,2).^2 );

% angular window (aka cos(\phi))
cwin = dot(xd,[cos(n0); sin(n0)]',2) ./ dr;

% amplitude correction w.r.t. reference line [Verhejen]
g0=sqrt( (xref(2)-x0(2,1)) / (xref(2) - xs(2)) );

% secondary source selection
if(lssel)
    idx = find(cwin>0);
else
    idx = 1:length(x0);
end

% calculate driving function only for active speakers
for n=idx
    D(n) = g0 .* sqrt(k/(2*pi*j)) * cwin(n) .* 1./dr(n).^(1/2) .* exp(-j*k*dr(n));
    %D(n) = g0 .* cwin(n) .* 1./dr(n).^(1/2) .* exp(-j*k*dr(n));
end