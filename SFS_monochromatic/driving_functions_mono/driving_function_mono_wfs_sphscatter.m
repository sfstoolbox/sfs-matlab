function D = driving_function_mono_wfs_sphscatter(x0,n0,Bl,f,xq,conf)

%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 6;
narginchk(nargmin,nargmax);
isargmatrix(x0,n0);
isargpositivescalar(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end
if nargin == nargmin
  xq = [0, 0, 0];
end
isargposition(xq);

%% ===== Configuration ==================================================
c = conf.c;
dimension = conf.dimension;
driving_functions = conf.driving_functions;
Nse = conf.scattering.Nse;
showprogress = conf.showprogress;

%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain


% apply shift to the center of expansion xq
x = x0(:,1)-xq(1);
y = x0(:,2)-xq(2);
z = x0(:,3)-xq(3);

% conversion to spherical coordinates
r = sqrt(x.^2 + y.^2 + z.^2);
phi = atan2(y,x);
theta = asin(z./r);

% frequency depended stuff
omega = 2*pi*f;
k = omega/c;
kr = k.*r;

% gradient in spherical coordinates
Gradr = zeros(size(x0,1),1);
Gradphi = zeros(size(x0,1),1);
Gradtheta = zeros(size(x0,1),1);

% directional weights for conversion spherical gradient into carthesian
% coordinates + point product with normal vector n0 (directional derivative
% in cartesian coordinates)
Sn0r     =  cos(theta).*cos(phi).*n0(:,1)...
         +  cos(theta).*sin(phi).*n0(:,2)...
         +  sin(theta)          .*n0(:,3);
Sn0phi   = -sin(phi)            .*n0(:,1)...
         +  cos(phi)            .*n0(:,2);
Sn0theta =  sin(theta).*cos(phi).*n0(:,1)...
         +  sin(theta).*sin(phi).*n0(:,2)...
         -  cos(theta)          .*n0(:,3);

% indexing the expansion coefficients
L = (Nse + 1).^2;
l = 0;
if strcmp('2.5D',dimension)
  if (strcmp('point_source',driving_functions))
    for n=0:Nse
      h_prime = k.*sphbesselh_derived(n,2,kr);
      h = sphbesselh(n, 2, kr);
      for m=-n:n
        l = l + 1;
        Ynm = sphharmonics(n,m, theta, phi);
        Gradr   = Gradr   +       ( Bl(l).*h_prime.*Ynm);
        Gradphi = Gradphi + 1./r.*( Bl(l).*h.*1j.*m.*Ynm );
        %Gradtheta = 0;
      end
      if showprogress, progress_bar(l,L); end % progress bar
    end
    % directional gradient + time reversion (conjugate complex)
    D = Sn0r.*Gradr + Sn0phi.*Gradphi + Sn0theta.*Gradtheta;
    D = -conj(D);
  else
    error(['%s: %s, this type of driving function is not implemented ', ...
      'for a 2.5D point source.'],upper(mfilename),driving_functions);
  end
else
  error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
