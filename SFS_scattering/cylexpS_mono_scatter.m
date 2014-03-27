function Bl = cylexpS_mono_scatter(Al, R, sigma, f, conf)

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
Nce = conf.scattering.Nce;

%% ===== Computation ====================================================
k = 2*pi*f/conf.c;
kR = k.*R;

if isinf(sigma)
  T = @(x) -besselj(x,kR)./besselh(x,2,kR);
elseif sigma == 0
  T = @(x) -besselj_derived(x,kR)./besselh_derived(x,2,kR);
else
  T = @(x) -(k.*besselj_derived(x,kR)+sigma.*besselj(x,kR)) ...
    ./(k.*besselh_derived(x,2,kR)+sigma.*besselh(x,2,kR));
end

L = 2*Nce+1;
Bl = zeros(L,1);
l = 0;
for n=-Nce:Nce
  l = l+1;
  Bl(l) = T(n).*Al(l);
  if showprogress, progress_bar(l,L); end % progress bar
end

end

