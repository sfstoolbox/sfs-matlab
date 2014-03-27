function [P] = sound_field_mono_sphbasis(Anm, Bn, Ynm,conf)

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
Nse = conf.scattering.Nse;

%% ===== Computation ====================================================
L = (Nse + 1).^2;
l = 0;

P = zeros(size(Bn{1}));
for n=0:Nse
  for m=-n:n
    l=l+1;    
    if showprogress, progress_bar(l,L); end  % progress bar
    P = P + Anm(l)*(Bn{n+1}.*Ynm{l});  
  end
end

end

