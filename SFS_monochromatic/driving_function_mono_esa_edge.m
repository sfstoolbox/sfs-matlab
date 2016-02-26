function D = driving_function_mono_esa_edge(x0,xs,src,f,conf)
% DRIVING_FUNCTION_MONO_ESA_EDGE 

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
isargsecondarysource(x0);
isargxs(xs);
isargpositivescalar(f);
isargchar(src);
isargstruct(conf);

%% ===== Configuration ==================================================
c = conf.c;
Xc = conf.secondary_sources.center;
alpha = conf.secondary_sources.alpha;

%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

if numel(alpha) ==2;
  phishift = alpha(1);
  alpha = alpha(2) - alpha(1);
else
  phishift = 0; 
end

% secondary source positions
x00 = bsxfun(@minus,x0(:,1:3),Xc);
[phi0, ~,r0] = cart2sph(x00(:,1),x00(:,2),x00(:,3));
phi0 = mod(phi0 - phishift, 2*pi);
% virtual source
[phis, ~,rs] = cart2sph(xs(1)-Xc(1),xs(2)-Xc(2),xs(3)-Xc(3));
phis = mod(phis - phishift, 2*pi);
% frequency depended stuff
omega = 2*pi*row_vector(f);  % [1 x Nf]
k = omega./c;  % [1 x Nf]
kr0 = r0 * k;  % [N0 x Nf]
krs = rs .* k;  % [1 x Nf]

Nce = 2*ceil(max(kr0));  % TODO: put this somewhere else

D = zeros(size(kr0));
epsilon = ones(Nce+1);
epsilon(1) = 2;
for n = 0:Nce
  nu = n*pi/alpha;
  factor = nu./epsilon(n+1) .* sin(nu*phis) * cos(nu*phi0) ./r0;

  select = r0 <= rs;
  
  D(select,:) = D(select,:) + factor(select) .* besselj(nu, kr0(select,:)) .*  ...
    besselh(nu, 2, krs);
  D(~select,:) = D(~select,:) + factor(~select) .* besselj(nu, krs) .* ...
    besselh(nu, 2, kr0(~select,:));
end

select = (-sin(phi0+phishift).*x0(:,4) + cos(phi0+phishift).*x0(:,5)) > 0;
D(select) = -D(select);
D = -j*D*pi/alpha;
