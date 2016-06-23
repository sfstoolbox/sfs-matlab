function [EF, EFm] = circexp_mono_translation(xt, mode, Nce, f, conf)
%CIRCEXP_MONO_TRANSLATION compute circular translation coefficients 
%(multipole re-expansion)
%
%   Usage: [EF, EFm] = circexp_mono_translation(xt, mode, Nce, f, conf)
%
%   Input parameters:
%       xt          - translatory shift [1x3] / m                    
%       mode        - 'RS' for regular-to-singular reexpansion
%                     'RR' for regular-to-regular reexpansion
%                     'SR' for singular-to-regular reexpansion
%                     'SS' for singular-to-singular reexpansion
%       f           - frequency / Hz
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       EF          - circular re-expansion coefficients for t
%       EFm         - circular re-expansion coefficients for -t
%  
%  CIRCEXP_MONO_TRANSLATION(t, mode, f, conf) computes the circular re-expansion
%  coefficients to perform as translatory shift of circular basis function.
%  Multipole Re-expansion computes the circular basis function for a shifted
%  coordinate system (x+t) based on the original basis functions for (x). 
%
%              \~~ inf
%  E (x + t) =  >         (E|F)   (xt) F (x)
%   n          /__ l=-inf      l,n      l
%
%  where {E,F} = {R,S}. R denotes the regular circular basis function, while
%  S symbolizes the singular circular basis function. Note that (S|S) and 
%  (S|R) are equivalent to (R|R) and (R|S), respectively.
%
%  see also: circexp_mono_ps, circexp_mono_pw
 
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
nargmin = 5;
nargmax = 5;
narginchk(nargmin,nargmax);
isargposition(xt);
isargchar(mode);
isargpositivescalar(Nce,f);

%% ===== Configuration ==================================================
c = conf.c;

%% ===== Variables ======================================================
% convert t into spherical coordinates
rt = norm(xt(1:2));
phit = atan2(xt(2),xt(1));

% frequency
k = 2*pi*f/c;
kr = k*rt;

L = 2*Nce+1;
EF = zeros(L,L);
EFm = EF;

% select suitable basis function
if strcmp('RR', mode) || strcmp('SS', mode)
  circbasis = @(nu) besselj(nu, kr);
elseif strcmp('SR', mode) || strcmp('RS', mode)
  circbasis = @(nu) besselh(nu, 2, kr);
else
  error('unknown mode:');
end

%% ===== Computation ====================================================
s = 0;
for n=-Nce:Nce
  s = s+1;
  l = 0;
  for m=-Nce:Nce
    l = l+1;
    EFm(s,l) = circbasis(n-m) .* exp(-1j.*(n-m).*phit);
    EF(s,l) = EFm(s,l) .* (-1)^(n-m);
  end
end

end