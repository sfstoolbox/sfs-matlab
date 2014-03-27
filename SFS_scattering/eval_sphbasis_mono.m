function [Jn, H2n, Ynm]  = eval_sphbasis_mono(r,theta,phi,k,conf)


%% ===== Configuration ==================================================
% Plotting result
Nse = conf.scattering.Nse;
showprogress = conf.showprogress;

%% ===== Computation ====================================================
kr = k.*r;  % argument of bessel functions

NJ = Nse + 1;
L = (NJ).^2;

Jn = cell(NJ,1);
H2n = cell(NJ,1);
Ynm = cell(L,1);

for n=0:Nse
  Jn{n+1} = sphbesselj(n,kr);
  H2n{n+1} = Jn{n+1} - 1j*sphbessely(n,kr);  
  for m=0:n
    l_plus = (n + 1).^2 - (n - m);
    l_minus = (n + 1).^2 - (n + m);
    if showprogress, progress_bar(l_plus,L); end  % progress bar
    % spherical harmonics (caution: symmetry relation depends on definition)
    Ynm{l_plus} = sphharmonics(n,m,theta,phi);
    Ynm{l_minus} = conj(Ynm{l_plus});
  end
end

end

