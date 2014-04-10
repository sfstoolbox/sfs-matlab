function [Alplus, SRplus, Alminus, SRminus] = cylexpSR_mono(Bl, xq, f, conf)
%Regular-To-Regular Cylindrical Reexpansion (Translatory shift of Expansion)
%
%   Usage: [Alplus, RRplus, Alminus, RRminus] = cylexpRR_mono(Al, xq, f, conf)
%
%   Input parameters:
%       Bl          - original singular cylindrical expansion coefficients [nx1]    
%       xq          - original expansion center coordinate
%       f           - frequency
%
%   Output parameters:
%       Alplus      - regular cylindrical expansion coefficients [nx1] 
%       SRplus      - singular-to-regular cylindrical reexpansion coefficients [nxn]
%       Alminus     - regular cylindrical expansion coefficients [nx1] 
%       SRminus     - singular-to-regular cylindrical reexpansion coefficients [nxn]     

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
isargposition(xq);
if nargin<nargmax
  conf = SFS_config;
else
  isargstruct(conf);
end

%% ===== Configuration ==================================================
showprogress = conf.showprogress;
Nce = conf.scattering.Nce;

%% ===== Computation ====================================================
% convert (xq) into spherical coordinates
r = sqrt(sum(xq(1:2).^2));
phi = atan2(xq(2),xq(1));

% frequency
k = 2*pi*f/conf.c;
kr = k*r;

L = 2*Nce+1;
SRplus = zeros(L,L);
SRminus = zeros(L,L);

s = 0;
for n=-Nce:Nce
  s = s+1;
  l = 0;
  for m=-Nce:Nce
    l = l+1;
    SRplus(s,l) = besselh(n-m,2,kr) .* exp(-1j.*(n-m).*phi);
    SRminus(s,l) = SRplus(s,l) .* -1.^(n-m);  
  end
  if showprogress, progress_bar(s,L); end  % progress bar
end

Alplus = SRplus*Bl;
Alminus = SRminus*Bl;

end

