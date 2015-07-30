function Anm = sphexp_mono_ps(xs, mode, Nse, f, xq, conf)
%Regular/Singular Spherical Expansion of Point Source
%
%   Usage: Anm = sphexp_mono_ps(xs, mode, Nse, f, xq, conf)
%
%   Input parameters:
%       xs          - position of point source
%       mode        - 'R' for regular, 'S' for singular
%       Nse         - maximum order of spherical basis functions
%       f           - frequency [m x 1] or [1 x m]
%       xq          - optional expansion center coordinate 
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       Anm         - regular Spherical Expansion Coefficients [(Nse+1)^2 x m]
%
%   SPHEXP_MONO_PS(xs, mode, f, Nse, xq, conf) computes the regular/singular 
%   Spherical Expansion Coefficients for a point source at xs. The expansion 
%   will be done around the expansion coordinate xq:
%
%   Regular Expansion:
%                \~~ oo  \~~   n   m  m
%   p    (x,f) =  >       >       A  R  (x-x ) 
%    ps,R        /__ n=0 /__ m=-n  n  n     q
%
%   with the expansion coefficients (Gumerov2004, eq. 3.2.2):
%    m               -m
%   A  = -i  . k  . S  (x  - x )
%    n               n   s    q
%
%   Singular Expansion:
%                \~~ oo  \~~   n   m  m
%   p    (x,f) =  >       >       B  S  (x-x ) 
%    ps,S        /__ n=0 /__ m=-n  n  n     q
%   
%   with the expansion coefficients (Gumerov2004, eq. 3.2.2):
%    m               -m
%   B  = -i  . k  . R  (x  - x )
%    n               n   s    q
%
%   The coefficients are stored in linear arrays with index l resulting from 
%   m and n:
% 
%         m         m               2
%   A  = A  ; B  = B  with l = (n+1)  - (n - m)
%    l    n    l    n
%
%   References:
%       Gumerov,Duraiswami (2004) - "Fast Multipole Methods for the 
%                                    Helmholtz Equation in three 
%                                    Dimensions", ELSEVIER
%
%   see also: sphexp_access sphexp_index sphbasis_mono

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************

%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 6;
narginchk(nargmin,nargmax);
isargposition(xs);
isargpositivescalar(Nse);
isargvector(f);
if nargin < nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end
if nargin == nargmin
  xq = [0,0,0];
else
  isargposition(xq);
end

%% ===== Configuration ==================================================
c = conf.c;

%% ===== Variables ======================================================
% convert (xs-xq) into spherical coordinates
r = sqrt(sum((xs-xq).^2));
phi = atan2(xs(2)-xq(2),xs(1)-xq(1));
theta = asin((xs(3)-xq(3))/r);

% frequency
k = 2.*pi.*row_vector(f)./c;
kr = k.*r;
Nf = length(kr);

% select suitable basis function
if strcmp('R', mode)
  sphbasis = @(nu) sphbesselh(nu,2,kr);
elseif strcmp('S', mode)
  sphbasis = @(nu) sphbesselj(nu, kr);
else
  error('unknown mode:');
end

%% ===== Computation ====================================================
Anm = zeros( (Nse + 1).^2 , Nf);
for n=0:Nse
  cn = -1j.*k.*sphbasis(n);
  for m=0:n    
    % spherical harmonics: conj(Y_n^m) = Y_n^-m (Gumerov2004, eq. 2.1.59)
    Ynm = sphharmonics(n,m, theta, phi);
    % -m
    v = sphexp_index(-m,n);  
    Anm(v,:) = cn.*Ynm;
    % +m
    v = sphexp_index(m,n); 
    Anm(v,:) = cn.*conj(Ynm);
  end
end

end

