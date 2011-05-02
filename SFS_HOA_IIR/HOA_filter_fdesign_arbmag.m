function [OUT, X ,F, Gd, F2] = HOA_filter_fdesign_arbmag(E, config, gd)


%  b = fdesign.arbmagnphase('n,f,h',config.N,config.f3, abs(E(:,14)));
%  Hd = design(b, 'freqsamp');
%  OUT = Hd;
%  plot(20*log10(abs(E(:,14))));
% end
%  

OUT = zeros(config.N + 1, config.nLS);

for n = 0:1:config.nLS-1

    m = E(:,n+1).'; % coefficeints of all LS
    b = fdesign.arbmag('n,f,a',config.N,config.f3, abs(m));
    validstructures(b)
    Hd = design(b,'freqsamp','filterstructure','dffir');
    OUT(:,n+1) =  Hd.Numerator;
    OUT(:,n+1) = circshift(OUT(:,n+1),[-round(sum(gd(:,n+1))) 0]);
end
 

 [X,F] = freqz( OUT(:,config.cLS), 1, config.n, config.fs); % frequency response of FIR Filter
 [Gd,F2] = grpdelay( OUT(:,config.cLS), 1, config.n, config.fs); % group delay of FIR Filter
 



