function [jn, h2n, Ynm] = sphbasis_mono(r,theta,phi,Nse,k,conf)
%Evaluate spherical basis functions for given input arguments
%
%   Usage: [jn, h2n, Ynm] = sphbasis_mono(r,theta,phi,Nse,k,conf)
%
%   Input parameters:
%       r           - distance from origin / m [n1 x n2 x ...]
%       theta       - elevation angle / rad [n1 x n2 x ...]
%       phi         - azimuth angle / rad [n1 x n2 x ...]
%       Nse         - maximum order of spherical basis functions
%       k           - wave number
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       jn          - cell array of spherical bessel functions
%       h2n         - cell array of spherical hankel functions of 2nd kind
%       Ynm         - cell array of spherical harmonics
%
%   SPHBASIS_MONO(r,theta,phi,k,conf) computes spherical basis functions for
%   the given arguments r, theta and phi. r, theta and phi can be of arbitrary
%   (but same) size. Output will be stored in cell arrays (one cell entry for 
%   each order) of length Nse+1 for jn and h2n. For Ynm the lenght is 
%   (Nse+1).^2. The coefficients of Ynm are stored with the linear index l 
%   resulting from the order m and the degree n of the spherical harmonics: 
%         m                 2
%   Y  = Y  ; with l = (n+1)  - (n - m)
%    l    n
%
%   see also: sphbasis_mono_grid sphharmonics sphbesselj sphbesselh

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
nargmin = 5;
nargmax = 6;
narginchk(nargmin,nargmax);
isargpositivescalar(Nse);
isargequalsize(r,phi,theta);
isargscalar(k);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end

%% ===== Configuration ==================================================
showprogress = conf.showprogress;

%% ===== Computation ====================================================
kr = k.*r;  % argument of bessel functions

NJ = Nse + 1;
L = (NJ).^2;

jn = cell(NJ,1);
h2n = cell(NJ,1);
Ynm = cell(L,1);

for n=0:Nse
  jn{n+1} = sphbesselj(n,kr);
  h2n{n+1} = jn{n+1} - 1j*sphbessely(n,kr);  
  for m=0:n
    l_plus = (n + 1).^2 - (n - m);
    l_minus = (n + 1).^2 - (n + m);
    if showprogress, progress_bar(l_plus,L); end  % progress bar
    % spherical harmonics (caution: symmetry relation depends on definition)
    Ynm{l_plus} = sphharmonics(n,m,theta,phi);
    Ynm{l_minus} = conj(Ynm{l_plus});
  end
end

end

