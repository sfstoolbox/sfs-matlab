function Bl = sphscatter_mono(Al, R, sigma, f, conf)

%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
isargvector(Al);
isargscalar(sigma);
isargpositivescalar(f,R);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end

%% ===== Configuration ==================================================
showprogress = conf.showprogress;
Nse = conf.scattering.Nse;

%% ===== Computation ====================================================
k = 2*pi*f/conf.c;
kR = k.*R;

L = (Nse + 1).^2;
Bl = zeros(L,1);

if isinf(sigma)
  T = @(x) -sphbesselj(x,kR)./sphbesselh(x,2,kR);
elseif sigma == 0
  T = @(x) -sphbesselj_derived(x,kR)./sphbesselh_derived(x,2,kR);
else
  T = @(x) -(k.*sphbesselj_derived(x,kR)+sigma.*sphbesselj(x,kR)) ...
    ./(k.*sphbesselh_derived(x,2,kR)+sigma.*sphbesselh(x,2,kR));
end

l = 0;
for n=0:Nse  
  fac = T(n);
  for m=-n:n
    l = l+1;
    % coefficients
    Bl(l) = fac.*Al(l);
  end
  if showprogress, progress_bar(l,L); end % progress bar
end

end

