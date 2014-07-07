function Al = cylexpR_mono_ls(xs,f,xq,conf)
%Regular Cylindrical Expansion of Line Source
%
%   Usage: Al = cylexpR_mono_pw(xs,f,xq,conf)
%
%   Input parameters:
%       xs          - position of line source 
%       f           - frequency
%       xq          - optional expansion center
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       Al          - regular cylindrical Expansion Coefficients
%
%   CYLEXPR_MONO_PW(nk,xq,f,conf) computes the regular cylindrical
%   expansion coefficients for a line source. The expansion will be done 
%   around the expansion coordinate xq:
%
%              \~~  oo       
%   p  (x,f) =  >        A  R  (x-x ) 
%    pw        /__ n=-oo  n  n     q
%
%   with the cylyndrical expansion coefficients:
%
%             n
%   A  = 4pi i  exp  (-i*n*phi  )
%    n                        pw
%
%   References:
%       Gumerov,Duraiswami (2004) - "Fast Multipole Methods for the 
%                                    Helmholtz Equation in three 
%                                    Dimensions", ELSEVIER
%
%   see also: eval_cylbasis_mono

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
  xq = [0, 0, 0];
else
  isargposition(xq);
end

%% ===== Configuration ==================================================
showprogress = conf.showprogress;
Nce = conf.scattering.Nce;
timereverse = conf.scattering.timereverse;
xref = conf.xref;
c = conf.c;

%% ===== Computation ====================================================
% convert (xs-xq) into cylindrical coordinates
r = sqrt((xs(1)-xq(1)).^2 + (xs(2)-xq(2)).^2);
phi = atan2(xs(2)-xq(2),xs(1)-xq(1));

% frequency
k = 2.*pi.*f./conf.c;
kr = k.*r;

if (timereverse)
  H = @(x) conj(1j/4*besselh(x,2,kr));
else
  H = @(x) 1j/4*besselh(x,2,kr);
end

L = 2*Nce+1;
Al = zeros(L,1);
l = 0;
for n=-Nce:Nce
  l = l+1;
  Al(l) = H(n).*exp(-1j*n*phi);
  if showprogress, progress_bar(l,L); end  % progress bar
end

end

