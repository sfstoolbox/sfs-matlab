function Dm = driving_function_mono_nfchoa_cht_circexp(Pm, f, conf)
%DRIVING_FUNCTION_MONO_NFCHOA_CHT_CIRCEXP computes the circular harmonics
%transform of nfchoa driving functions for a sound field expressed by regular
%circular expansion coefficients.
%
%   Usage: D = driving_function_mono_nfchoa_cht_circexp(Pm, f, conf)
%
%   Input parameters:
%       Pm          - regular circular expansion coefficients of virtual
%                     sound field [n x m]
%       f           - frequency in Hz [m x 1] or [1 x m]
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       Dm          - circular harmonics transform of driving function signal
%                     [n x m]
%
%   DRIVING_FUNCTION_MONO_NFCHOA_CHT_CIRCEXP(Pm, f, conf) returns circular
%   harmonics transform of the NFCHOA driving function for a virtual sound
%   expressed by regular circular expansion coefficients and the frequency f.

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
isargmatrix(Pm);
isargvector(f);
isargstruct(conf);
if size(Pm,2) ~= length(f)
  error( '%s:number of rows in %s have to match length of %s', ...
    upper(mfilename), inputname(1), inputname(2) );
end

%% ===== Configuration ==================================================
c = conf.c;
R0 = conf.secondary_sources.size / 2;

%% ===== Variables ======================================================
Nce = (size(Pm, 1)-1)/2;

% frequency depended stuff
omega = 2*pi*row_vector(f);  % [1 x Nf]
k = omega./c;  % [1 x Nf]
kR0 = k.*R0;  % [1 x Nf]

%% ===== Computation ====================================================
% Calculate the circular harmonics of driving function
%
% Regular circular expansion of the sound field:
%          \~~ oo                 
% P(x,w) =  >        P  J (kr) . e^(+j m phi)
%          /__ m=-oo  m  m
%
% and 3D free field Green's Function:
%            \~~ oo                 
% G  (x,w) =  >        G  J (kr) . e^(+j m phi)
%  ls        /__ m=-oo  m  m
%
% with the regular expansion coefficients of Green's Function
% (see Gumerov2004, eq. 3.2.2):
%    m     -i      (2)
%   G  =  ----- . H   (kr0)
%    n      4      m

%                P
%        1        m
% D  = ------  --------
%  m   2pi r0    G
%                 m

Dm = zeros(size(Pm));
l = 0;
for m=-Nce:Nce
  l = l+1;
  Dm(l,:) = Pm(l,:) ./ besselh(m,2,kR0);
end
% factor from expansion of 2D free field Green's Function
Dm = Dm./(-1j/4);
% normalization due to size of circular array
Dm = Dm./2*pi*R0;
