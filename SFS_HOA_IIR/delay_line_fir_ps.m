function delay = delay_line_fir_ps(config)

delay=zeros(1,length(config.x0));
%amp=zeros(1,length(x0));

xd = bsxfun(@minus,config.x0',[config.xps(1) config.xps(2)]);
dr=sqrt( xd(:,1).^2 + xd(:,2).^2 );

%cwin = dot(xd,[cos(n0); sin(n0)]',2);

%g0=1;
%idx = find(cwin>0);


for n = 0:length(config.x0)-1   
    %amp(n) = g0 * cwin(n) .* 1./dr(n).^(3/2);
    delay(n+1) = 1/config.c * dr(n+1);
end

%delay = delay + max(abs(delay)); 
delay = round(delay.*config.fs);

end