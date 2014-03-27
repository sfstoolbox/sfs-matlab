function [Jn, H2n, Yn]  = eval_cylbasis_mono(r,phi,k,conf)


%% ===== Configuration ==================================================
% Plotting result
Nce = conf.scattering.Nce;
showprogress = conf.showprogress;

%% ===== Computation ====================================================
kr = k.*r;  % argument of bessel functions

L = 2*Nce + 1;
Jn = cell(L,1);
H2n = cell(L,1);
Yn = cell(L,1);

l = 0;
for n=-Nce:Nce
  l = l + 1;
  Jn{l} = besselj(n,kr);
  H2n{l} = Jn{l} - 1j*bessely(n,kr);
  Yn{l} = exp(1j*n*phi);
  if showprogress, progress_bar(l,L); end  % progress bar
end

end

