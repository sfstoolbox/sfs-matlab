function [OUT, X ,F, Gd, F2] = HOA_filter_fdesign_arbmagnphase(E, config)

OUT = zeros(config.N + 1, config.nLS);

for n = 0:1:config.nLS-1

   m = E(:,n+1).';
   b = fdesign.arbmagnphase('n,f,h',config.N,config.f3,abs(m));
   validstructures(b)
   Hd = design(b, 'freqsamp','filterstructure','dffir');
   OUT(:,n+1) = Hd.Numerator;
   
end

 [X,F] = freqz( OUT(:,config.cLS), 1, config.n, config.fs); % frequency response of FIR Filter
 [Gd,F2] = grpdelay( OUT(:,config.cLS), 1, config.n, config.fs); % group delay of FIR Filter
 



