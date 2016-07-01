function Dnm = driving_function_mono_nfchoa_sht_sphexp(Pnm, f, conf)
%DRIVING_FUNCTION_MONO_NFCHOA_SHT_SPHEXP computes the spherical harmonics
%transform of nfchoa driving functions for a sound field expressed by regular
%spherical expansion coefficients.
%
%   Usage: D = driving_function_mono_nfchoa_sht_sphexp(Pnm, f, conf)
%
%   Input parameters:
%       Pnm         - regular spherical expansion coefficients of virtual
%                     sound field [nxm]
%       f           - frequency in Hz [m x 1] or [1 x m]
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       Dnm         - regular spherical harmonics transform of driving
%                     function signal [n x m]
%
%   DRIVING_FUNCTION_MONO_NFCHOA_SHT_SPHEXP(Pnm, f, conf) returns spherical
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
R0 = conf.secondary_sources.size / 2;

%% ===== Variables ======================================================
Nse = sqrt(size(Pnm,1))-1;

% frequency depended stuff
omega = 2*pi*row_vector(f);  % [1 x Nf]
k = omega./c;  % [1 x Nf]
kR0 = k.*R0;  % [1 x Nf]

%% ===== Computation ====================================================
% Calculate the spherical harmonics of driving function
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

%             m
%            P
%  m          n
% D = ---------------------
%  n             (2)
%     -j k r0^2 h  (k r0)
%                n

Dnm = zeros(size(Pnm));
for n=0:Nse
  v = sphexp_index(-n:n, n); % m=-n:n
  Dnm(v,:) = bsxfun(@rdivide, Pnm(v,:), sphbesselh(n,2,kR0));
end
% order independent factor
Dnm = bsxfun(@rdivide, Dnm, -1j.*(R0.^2)*k);
