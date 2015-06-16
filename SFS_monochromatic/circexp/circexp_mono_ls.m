function ABm = circexp_mono_ls(xs, mode, Nce, f, xq, conf)
%regular/singular circular expansion  of line source
%
%   Usage: ABm = circexp_mono_ls(xs, mode, Nce, f, xq, conf)
%
%   Input parameters:
%       xs          - position of line source
%       mode        - 'R' for regular, 'S' for singular
%       Nce         - maximum order of circular basis functions
%       f           - frequency / Hz [1 x m] or [m x 1]
%       xq          - optional expansion center / m [1 x 3]
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       ABm         - regular cylindrical Expansion Coefficients
%
%   CIRCEXP_MONO_LS(xs, mode, Nce, f, xq, conf) computes the regular circular
%   expansion coefficients for a line source. The expansion will be done 
%   around the expansion coordinate xq:
%
%   Regular Expansion:
%                \~~ oo        
%   p    (x,f) =  >         A  R  (x-x ) 
%    ls,R        /__ n=-oo   n  n     q
%
%   with the expansion coefficients:
%    m   -i           
%   A  = ---  S (x -x ) 
%    n    4    n  s  q  
%
%   Singular Expansion:
%                \~~ oo      m  
%   p    (x,f) =  >         B  S  (x-x ) 
%    ls,R        /__ n=-oo   n  n     q
%
%   with the expansion coefficients):
%    m   -i           
%   A  = ---  R (x -x ) 
%    n    4    n  s  q

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
nargmin = 4;
nargmax = 6;
narginchk(nargmin,nargmax);
isargposition(xs);
isargchar(mode);
isargpositivescalar(Nce);
isargvector(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end
if nargin == nargmin
    xq = [0,0,0];
else
  isargposition(xq);
end

%% ===== Configuration ==================================================
c = conf.c;

%% ===== Variables ======================================================
% convert (xs-xq) into cylindrical coordinates
r = sqrt((xs(1)-xq(1)).^2 + (xs(2)-xq(2)).^2);
phi = atan2(xs(2)-xq(2),xs(1)-xq(1));

% frequency
k = 2.*pi.*row_vector(f)./c;
kr = k.*r;

% select suitable basis function
if strcmp('R', mode)
  circbasis = @(nu,z) besselh(nu,2,z);
elseif strcmp('S', mode)
  circbasis = @besselj;
else
  error('unknown mode:');
end

%% ===== Computation ====================================================
L = 2*Nce+1;
Nf = length(kr);
ABm = zeros(L,Nf);
l = 0;
for n=-Nce:Nce
  l = l+1;
  ABm(l,:) = -1j/4*circbasis(n,kr).*exp(-1j*n*phi);
end

end

