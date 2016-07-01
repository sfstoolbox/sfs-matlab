function Dm = driving_function_mono_nfchoa_cht_sphexp(Pnm, f, conf)
%DRIVING_FUNCTION_MONO_NFCHOA_CHT_SPHEXP computes the circular harmonics
%transform of nfchoa driving functions for a sound field expressed by regular
%spherical expansion coefficients.
%
%   Usage: D = driving_function_mono_nfchoa_cht_sphexp(Pnm, f, conf)
%
%   Input parameters:
%       Pnm         - regular spherical expansion coefficients of virtual
%                     sound field [nxm]
%       f           - frequency in Hz [m x 1] or [1 x m]
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       Dm          - regular circular harmonics transform of driving
%                     function signal
%
%   DRIVING_FUNCTION_MONO_NFCHOA_CHT_SPHEXP(Pnm, f, conf) returns circular
%   harmonics transform of the NFCHOA driving function for a virtual sound
%   expressed by regular spherical expansion coefficients and the frequency f.

%*****************************************************************************
% Copyright (c) 2010-2016 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2016 Institut fuer Nachrichtentechnik                   *
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
nargmin = 3;
nargmax = 3;
narginchk(nargmin,nargmax);
isargmatrix(Pnm);
isargsquaredinteger(size(Pnm,1));
isargvector(f);
isargstruct(conf);
if size(Pnm,2) ~= length(f)
  error( '%s:number of rows in %s have to match length of %s', ...
    upper(mfilename), inputname(1), inputname(2) );
end

%% ===== Configuration ==================================================
c = conf.c;
Xc = conf.secondary_sources.center;
R0 = conf.secondary_sources.size / 2;
xref = conf.xref - Xc;
rref = norm( xref );  % reference radius
thetaref = asin( xref(3)./rref);  % reference elevation angle

%% ===== Variables ======================================================
Nse = sqrt(size(Pnm,1))-1;

% frequency depended stuff
omega = 2*pi*row_vector(f);  % [1 x Nf]
k = omega./c;  % [1 x Nf]
kR0 = k.*R0;  % [1 x Nf]
krref = k.*rref;  % [1 x Nf]

%% ===== Computation ====================================================
% Calculate the circular harmonics of driving function
%
% Regular spherical expansion of the sound field:
%          \~~   N \~~   n   m           m
% P(x,w) =  >       >       P  j (kr) . Y  (theta, phi)
%          /__ n=0 /__ m=-n  n  n        n 
%
% and 3D free field Green's Function:
%             \~~ oo  \~~   n   m             m
% G  (x0,f) =  >       >       G  . j (kr) . Y  (theta, phi)
%  ps         /__ n=0 /__ m=-n  n    n        n
%
% with the regular expansion coefficients of Green's Function
% (see Gumerov2004, eq. 3.2.2):
%    m               (2)       -m
%   G  = -i  . k  . h   (kr0) Y  (pi/2, 0)
%    n               n         n

% initialize empty driving signal
Dm = zeros(2*Nse+1, size(Pnm,2));
l = 0;
if (rref == 0)
  % --- Xref == Xc -------------------------------------------------
  %                 m
  %                P
  %        1        |m|
  % D  = ------  --------
  %  m   2pi r0     m
  %                G
  %                 |m|
  for m=-Nse:Nse
    l = l+1;
    v = sphexp_index(m);  % n = abs(m)

    Dm(l,:) = Pnm(v,:) ./ ...
      ( -1i.*k.* sphbesselh(abs(m),2,kR0).*sphharmonics(abs(m),-m, 0, 0) );
  end
else
  % --- Xref ~= Xc --------------------------------------------------
  %
  % if the reference position is not in the middle of the array,
  % things get a 'little' more complicated
  %                __
  %               \     m             m
  %               /__  P  j (k rref) Y (thetaref, 0)
  %  n     1      l=|m| l  l          l
  % D = ------- ------------------------------------
  %  m  2pi r0     __
  %               \     m             m
  %               /__  G  j (k rref) Y (thetaref, 0)
  %               l=|m| n  l          l
  %

  % save some intermediate results   
  hn = zeros(size(Pnm));
  jn = hn;
  for n=0:Nse
    hn(n+1,:) = sphbesselh(n,2,kR0);  % [1 x Nf]
    jn(n+1,:) = sphbesselj(n,krref);  % [1 x Nf]
  end

  for m=-Nse:Nse
    Pm = zeros(1,length(f));
    Gm = Pm;
    for n=abs(m):Nse
      v = sphexp_index(m,n);
      % 
      factor = jn(n+1,:) .* sphharmonics(n, -m, thetaref, 0);
      % numerator
      Pm = Pm + Pnm(v,:) .* factor;
      % denominator
      Gm = Gm + (-1i*k) .* hn(n+1,:) .* sphharmonics(n, -m, 0, 0) .* factor;
    end
    l = l+1;
    Dm(v,:) = Pm./Gm;
  end
end
% normalization due to size of circular array
Dm = Dm./(2*pi*R0);