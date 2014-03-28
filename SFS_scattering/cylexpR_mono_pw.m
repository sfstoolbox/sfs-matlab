function Al = cylexpR_mono_pw(nk,f,xq,conf)
%Regular Cylindrical Expansion of Plane Wave
%
%   Usage: Al = cylexpR_mono_pw(nk,f,xq,conf)
%
%   Input parameters:
%       nk          - propagation direction of plane wave 
%       f           - frequency
%       xq          - optional expansion center
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       Al          - regular cylindrical Expansion Coefficients
%
%   CYLEXPR_MONO_PW(nk,xq,f,conf) computes the regular cylindrical
%   expansion coefficients for a plane wave. The expansion will be done a
%   round the expansion coordinate xq:
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
isargposition(nk);
isargpositivescalar(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end
if nargin == nargmin
  xq = [0, 0, 0];
end
isargposition(xq);

%% ===== Configuration ==================================================
showprogress = conf.showprogress;
Nce = conf.scattering.Nce;
timereverse = conf.scattering.timereverse;
xref = conf.xref;
c = conf.c;

%% ===== Computation ====================================================
if (timereverse)
  n_sign = 1;
else
  n_sign = -1;
end

% convert nk into cylindrical coordinates
phi = atan2(nk(2),nk(1));

% delay of plane wave to reference point
nk = nk./vector_norm(nk,2);
delay = 2*pi*f/c*(nk*(xq-xref)');

L = 2*Nce+1;
Al = zeros(L,1);
l = 0;
for n=-Nce:Nce
  l = l+1;
  Al(l) = (1j)^(n_sign*n).*exp(-1j*n*phi).*exp(-1j*delay);
  if showprogress, progress_bar(l,L); end  % progress bar
end

end

