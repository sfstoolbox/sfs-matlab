function [Jn, H2n, Yn]  = cylbasis_mono(r,phi,k,conf)
%Evaluate cylindrical basis functions for given input arguments
%
%   Usage: [Jn, H2n, Yn]  = cylbasis_mono(r,phi,k,conf)
%
%   Input parameters:
%       r           - distance from z-axis in cylindrical coordinates 
%       phi         - azimuth angle in cylindrical coordinates
%       k           - wave number
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       Jn          - cell array of cylindrical bessel functions
%       H2n         - cell array of cylindrical hankel functions of 2nd kind
%       Yn          - cell array of cylindrical harmonics
%
%   CYLBASIS_MONO(r,phi,k,conf) computes cylindrical basis functions for
%   the given arguments r and phi. r and phi can be of arbitrary (but same)
%   size. Output will be stored in cell arrays (one cell entry for each order)
%   of length 2*conf.scattering.Nce+1 . Each cell array entry contains a 
%   matrix of the same size as r and phi.
%
%   References:
%       Williams (1999) - "Fourier Acoustics", ACADEMIC PRESS
%
%   see also: cylbasis_mono_XYZgrid besselj besselh

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
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargequalsize(r,phi);
isargscalar(k);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end

%% ===== Configuration ==================================================
% Plotting result
Nce = conf.scattering.Nce;
showprogress = conf.showprogress;

%% ===== Computation ====================================================
kr = k.*r;  % argument of bessel functions

L = 2*Nce + 1;
Jn = cell(L,1);
H2n = cell(L,1);
Yn = cell(L,1);

l = 0;
for n=-Nce:0
  % negative n
  l = l + 1;
  Jn{l} = besselj(n,kr);
  H2n{l} = Jn{l} - 1j*bessely(n,kr);
  Yn{l} = exp(1j*n*phi);
  % positive n
  Jn{L-l+1} = (-1)^n*Jn{l};
  H2n{L-l+1} = (-1)^n*H2n{l};
  Yn{L-l+1} = conj(Yn{l});  
  if showprogress, progress_bar(n+Nce,Nce); end  % progress bar
end

end

