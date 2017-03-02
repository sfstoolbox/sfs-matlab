function Bm = circexp_mono_scatter(Am,R,sigma,f,conf)
%CIRCEXP_MONO_SCATTER computes the singular circular expansion coefficients of
%cylinder-scattered field
%
%   Usage: Bm = circexp_mono_scatter(Al,R,sigma,f,conf)
%
%   Input parameters:
%       Am          - regular circular expansion of incident field [n x Nf]
%       R           - radius of cylinder / m
%       sigma       - complex admittance of scatterer
%       f           - frequency / Hz [1 x Nf] or [Nf x 1]
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       Bm          - singular circular expansion coefficients of
%                     scattered field [n x Nf]
%
%   CIRCEXP_MONO_SCATTER(Am,R,sigma,f,conf) computes the singular
%   cylindrical expansion coefficients of a field resulting from a scattering
%   of an incident field at a cylindrical. Incident field is descriped by
%   regular expansion coefficients (expansion center is expected to be at the
%   center of the cylinder xq):
%
%               \~~ oo
%   p   (x,f) =  >      A  R  (x-x )
%    ind        /__ n=0  n  n     q
%
%   The scattered field is descriped by singular expansion coefficients,
%   expanded around the center of the cylinder x0.
%
%               \~~ oo
%   p   (x,f) =  >      B  S  (x-x )
%    sca        /__ n=0  n  n     q
%
%   Due to the boundary conditions on the surface of the cylinder the
%   coefficients are related by:
%
%          k  . J' (kR) + sigma  . H (kR)
%                n                 n
%   B  = - ------------------------------ . A
%    n     k  . H' (kR) + sigma  . H (kR)    n
%                 n                 n
%
%   where k = 2*pi*f/c.
%
%   see also: circexp_mono_pw circexp_mono_ls

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
isargmatrix(Am);
isargscalar(sigma);
isargpositivescalar(R);
isargvector(f);

%% ===== Configuration ==================================================
c = conf.c;

%% ===== Variables ======================================================
% frequency
k = 2.*pi.*row_vector(f)./c;
kR = k.*R;
% select suitable transformation function
if isinf(sigma)
    T = @(x) -besselj(x,kR)./besselh(x,2,kR);
elseif sigma == 0
    T = @(x) -besselj_derived(x,kR)./besselh_derived(x,2,kR);
else
    T = @(x) -(k.*besselj_derived(x,kR)+sigma.*besselj(x,kR)) ...
               ./(k.*besselh_derived(x,2,kR)+sigma.*besselh(x,2,kR));
end

%% ===== Computation ====================================================
Nce = (size(Am,1)-1)/2;
Bm = zeros(size(Am));
l = 0;
for n=-Nce:Nce
    l = l+1;
    Bm(l,:) = T(n).*Am(l,:);
end

end
