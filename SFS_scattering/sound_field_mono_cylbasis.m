function [P] = sound_field_mono_cylbasis(Anm, Bn, Yn,conf)

%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end

%% ===== Configuration ==================================================
% Plotting result
showprogress = conf.showprogress;
Nce = conf.scattering.Nce;

%% ===== Computation ====================================================
L = 2*Nce + 1;
l = 0;

P = zeros(size(Bn{1}));
for n=-Nce:Nce
  l=l+1;
  P = P + Anm(l)*(Bn{l}.*Yn{l});
  if showprogress, progress_bar(l,L); end  % progress bar
end

end

