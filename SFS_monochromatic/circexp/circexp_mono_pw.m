function Am = circexp_mono_pw(nk, Nce, f, xq, conf)
%regular circular expansion of plane wave
%
%   Usage: Am = circexp_mono_pw(nk, Nce, f, xq, conf)
%
%   Input parameters:
%       nk          - propagation direction of plane wave 
%       Nce         - maximum order of circular basis functions
%       f           - frequency / Hz [1 x m] or [m x 1]
%       xq          - optional expansion center
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       Am          - regular cylindrical Expansion Coefficients
%
%   CIRCEXP_MONO_PW(nk, Nce, f, xq, conf) computes the regular circular
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
nargmin = 3;
nargmax = 5;
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
c = conf.c;

%% ===== Variables ======================================================
% convert nk into cylindrical coordinates
phi = atan2(nk(2),nk(1));

% delay of plane wave to reference point
nk = nk./vector_norm(nk,2);
delay = 2*pi*row_vector(f)/c*(nk*xq.');

%% ===== Computation ====================================================
L = 2*Nce+1;
Nf = length(delay);
Am = zeros(L,Nf);
l = 0;
for n=-Nce:Nce
  l = l+1;
  Am(l,:) = (1j)^(-n).*exp(-1j*n*phi).*exp(-1j*delay);
end

end

