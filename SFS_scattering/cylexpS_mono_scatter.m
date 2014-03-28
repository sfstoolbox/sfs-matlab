function Bl = cylexpS_mono_scatter(Al, R, sigma, f, conf)
%Singular Spherical Expansion of cylinder-scattered field
%
%   Usage: Bl = cylexpS_mono_scatter(Al, R, sigma, f, conf)
%
%   Input parameters:
%       Al          - regular cylindrical expansion of incident field (cylexpR_*)                     
%       R           - radius of cylinder
%       sigma       - complex admittance of scatterer
%       f           - frequency in Hz
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       Bl          - singular cylindrical expansion coefficients of
%                     scattered field
%
%   CYLEXPS_MONO_SCATTER(Al, R, sigma, f, conf) computes the singular 
%   cylindrical expansion coefficients of a field resulting from a scattering 
%   of an incident field at a cylindrical. Incident field is descriped by 
%   regular expansion coefficients (expansion center is expected to be at the 
%   center of the cylinder xq):
%
%               \~~ oo     
%   p   (x,f) =  >      A  R  (x-x ) 
%    ind        /__ n=0  n  n     q
%
%   The scattered field is descriped by singular expansion coefficients,
%   expanded around the center of the cylinder x0. 
%
%               \~~ oo    
%   p   (x,f) =  >      B  S  (x-x ) 
%    sca        /__ n=0  n  n     q
%
%   Due to the boundary conditions on the surface of the cylinder the
%   coefficients are related by:
%
%          k  . J' (kR) + sigma  . H (kR)
%                n                 n        
%   B  = - ------------------------------ . A
%    n     k  . H' (kR) + sigma  . H (kR)    n
%                 n                 n
%
%   where k = 2*pi*f/c.
%
%   References:
%       Gumerov,Duraiswami (2004) - "Fast Multipole Methods for the 
%                                    Helmholtz Equation in three 
%                                    Dimensions", ELSEVIER
%
%   see also: cylexpR_mono_pw



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

