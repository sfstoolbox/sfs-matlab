function Al = sphexpR_mono_ps(xs,f,xq,conf)
%Regular Spherical Expansion of Plane Wave
%
%   Usage: Al = sphexpR_mono_ps(nk,f,x0,conf)
%
%   Input parameters:
%       xs          - position of 
%       f           - frequency
%       xq          - optional expansion coordinate 
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       Al          - regular Spherical Expansion Coefficients
%
%   SPHEXPR_MONO_PS(xs,f,xq,conf) computes the regular Spherical Expansion
%   Coefficients for a plane wave. The expansion will be done around the
%   expansion coordinate x0.  
%
%   References:
%       Gumerov,Duraiswami (2004) - "Fast Multipole Methods for the 
%                                    Helmholtz Equation in three 
%                                    Dimensions", ELSEVIER
%
%   see also:

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
nargmin = 2;
nargmax = 4;
narginchk(nargmin,nargmax);
isargposition(xs);
isargpositivescalar(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end
if nargin == nargmin
    xq = [0,0,0];
end
isargposition(xq);

%% ===== Configuration ==================================================
showprogress = conf.showprogress;
Nse = conf.scattering.Nse;
timereverse = conf.scattering.timereverse;

%% ===== Computation ====================================================
% convert (xs-xq) into spherical coordinates
r = sqrt((xs(1)-xq(1)).^2 + (xs(2)-xq(2)).^2 + (xs(3)-xq(3)).^2);
phi = atan2(xs(2)-xq(2),xs(1)-xq(1));
theta = asin((xs(3)-xq(3))/r);

% frequency
k = 2.*pi.*f./conf.c;
kr = k.*r;

if (timereverse)
  H = @(x) conj(1j*k*sphbesselh(x,2,kr));
else
  H = @(x) 1j*k*sphbesselh(x,2,kr);
end

L = (Nse + 1).^2;
Al = zeros(L,1);
for n=0:Nse
  Hn = H(n);
  for m=0:n    
    l_plus = (n + 1).^2 - (n - m);
    l_minus = (n + 1).^2 - (n + m);    
    % caution: symmetry relation depends on definition of spherical harmonics
    Ynm = sphharmonics(n,-m, theta, phi);  % spherical harmonics
    Al(l_plus) = Hn.*Ynm;
    Al(l_minus) = Hn.*conj(Ynm);
  end
  if showprogress, progress_bar(l_plus,L); end % progress bar
end

end

