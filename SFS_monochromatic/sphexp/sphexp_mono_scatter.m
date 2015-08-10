function Bnm = sphexp_mono_scatter(Anm, R, sigma, f, conf)
%SPHEXP_MONO_SCATTER compute the singular spherical expansion of a sphere-
%scattered sound field
%
%   Usage: Bnm = sphexp_mono_scatter(Anm, R, sigma, f, conf)
%
%   Input parameters:
%       Anm         - regular spherical expansion of incident sound field
%       R           - radius of sphere / m
%       sigma       - admittance of sphere (from 0(sound soft) to inf(sound hard))
%       f           - frequency / Hz [Nf x 1] or [1 x Nf]
%
%   Output parameters:
%       Bnm         - singular spherical expansion coefficients of
%                     scattered field
%
%   SPHEXP_MONO_SCATTER(Anm, R, sigma, f, conf) computes the singular spherical
%   expansion coefficients of a field resulting from a scattering of an incident
%   field at a sphere. Incident field is descriped by regular expansion
%   coefficients (expansion center is expected to be at the center of the sphere
%   xq):
%
%               \~~ oo  \~~   n   m  m
%   P   (x,f) =  >       >       A  R  (x-x )
%    ind        /__ n=0 /__ m=-n  n  n     q
%
%   The scattered field is descriped by singular expansion coefficients,
%   expanded around the center of the sphere xq.
%
%               \~~ oo  \~~   n   m  m
%   P   (x,f) =  >       >       B  S  (x-x )
%    sca        /__ n=0 /__ m=-n  n  n     q
%
%   Due to the boundary conditions on the surface of the sphere the
%   coefficients are related by:
%
%          k  . j' (kR) + sigma  . j (kR)
%    m            n                 n        m
%   B  = - ------------------------------ . A
%    n     k  . h' (kR) + sigma  . h (kR)    n
%                 n                 n
%
%   where k = 2*pi*f/c. The coefficients are stored in linear arrays with
%   index l resulting from m and n:
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
%   see also: sphexp_mono_ps, sphexp_mono_pw, sphexp_mono_timereverse

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
nargmax = 5;
narginchk(nargmin,nargmax);
isargvector(Anm);
L = length(Anm);
isargsquaredinteger(L);
isargscalar(sigma);
isargpositivescalar(R);
isargvector(f);
if nargin<nargmax
  conf = SFS_config;
else
  isargstruct(conf);
end

%% ===== Configuration ==================================================
c = conf.c;

%% ===== Variables ======================================================
% frequency
k = 2.*pi.*row_vector(f)./c;
kR = k.*R;

% select suitable transformation function
if isinf(sigma)
  T = @(x) -sphbesselj(x,kR)./sphbesselh(x,2,kR);
elseif sigma == 0
  T = @(x) -sphbesselj_derived(x,kR)./sphbesselh_derived(x,2,kR);
else
  T = @(x) -(k.*sphbesselj_derived(x,kR)+sigma.*sphbesselj(x,kR)) ...
    ./(k.*sphbesselh_derived(x,2,kR)+sigma.*sphbesselh(x,2,kR));
end

%% ===== Computation ====================================================
Bnm = zeros(size(Anm));
for n=0:(sqrt(L) - 1)
  v = sphexp_index(-n:n,n);
  Bnm(v,:) = T(n).*Anm(v,:);
end

end

