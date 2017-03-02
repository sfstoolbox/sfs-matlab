function D = driving_function_mono_nfchoa_sphexp(x0, Pnm,f,conf)
%computes the nfchoa driving functions for a sound field expressed by regular
%spherical expansion coefficients.
%
%   Usage: D = driving_function_mono_nfchoa_sphexp(x0,Pnm,f,conf)
%
%   Input parameters:
%       x0          - position of the secondary sources / m [N0x3]
%       Pnm         - regular spherical expansion coefficients of sound field
%                     [N x Nf]
%       f           - frequency / Hz [Nf x 1] or [1 x Nf]
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_NFCHOA_SPHEXP(x0,Pnm,f,conf) returns NFCHOA
%   driving signals for the given secondary sources, the virtual sound expressed
%   by regular spherical expansion coefficients and the frequency f.
%
%   see also: driving_function_mono_wfs_sphexp

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
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargmatrix(Pnm,x0);
isargvector(f);
isargstruct(conf);
if mod(sqrt(size(Pnm, 1)),1) ~= 0
  error(['%s: number of columns of Pnm (%s) is not the square of an', ...
    'integer.'], upper(mfilename), sqrt(size(Pnm, 1)));
end

%% ===== Configuration ==================================================
c = conf.c;
dimension = conf.dimension;
Xc = conf.secondary_sources.center;
R0 = conf.secondary_sources.size/2;

%% ===== Variables ======================================================
Nse = sqrt(size(Pnm, 1))-1;

% secondary source positions
x00 = bsxfun(@minus,x0(:,1:3),Xc);
[phi0, theta0, r0] = cart2sph(x00(:,1),x00(:,2),x00(:,3));

% frequency depended stuff
omega = 2*pi*row_vector(f);  % [1 x Nf]
k = omega./c;  % [1 x Nf]
kr0 = r0 * k;  % [N0 x Nf]

% initialize empty driving signal
D = zeros(size(kr0));  % [N0 x Nf]

%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

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
% with the regular expansion coefficients of reen's Function
% (see Gumerov2004, eq. 3.2.2):
%    m               (2)       -m
%   G  = -i  . k  . h   (kr0) Y  (pi/2, 0)
%    n               n         n

if strcmp('2D',dimension)
  % === 2-Dimensional ==================================================
  % Convert spherical Expansion to circular expansion (approximation)
  Pm = circexp_convert_sphexp(Pnm);
  % Compute driving function for circular expansion
  D = driving_function_mono_nfchoa_circexp(x0,Pm,f,conf);
  %
elseif strcmp('2.5D',dimension)
  % === 2.5-Dimensional ================================================
  % Spherical Harmonics of Driving Function
  Dm = driving_function_mono_nfchoa_cht_sphexp(Pnm,f,conf);
  
  l = 0;
  for m=-Nse:Nse
    l=l+1;
    % second line is for compensating possible difference between r0 and R0
    D = D + bsxfun(@times, Dm(l,:), exp(1i.*m.*phi0)).* ...
      bsxfun(@rdivide, sphbesselh(abs(m),2,kr0), sphbesselh(abs(m),2,k*R0));
  end
  %
elseif strcmp('3D',dimension)
  % === 3-Dimensional ==================================================
  % Spherical Harmonics of Driving Function
  Dnm = driving_function_mono_nfchoa_sht_sphexp(Pnm,f,conf);
  
  % Inverse Spherical Harmonics Transform
  for n=0:Nse
    Dn = 0;
    for m=-n:n
      v = sphexp_index(m, n);
      Dn = Dn + bsxfun(@times, Dnm(v,:), sphharmonics(n,m,theta0,phi0));
    end
    % compensating possible difference between r0 and R0
    D = D + Dn .* bsxfun(@rdivide, sphbesselh(n,2,kr0), sphbesselh(n,2,k*R0));
  end
  %
else
  error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
