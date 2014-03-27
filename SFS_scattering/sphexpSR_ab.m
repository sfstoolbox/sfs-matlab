function [a, b] = sphexpSR_ab(Nse,conf)
%
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
nargmin = 0;
nargmax = 2;
narginchk(nargmin,nargmax);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end
if nargin == nargmin
  Nse = conf.scattering.Nse;
end

%% ===== Configuration ==================================================
showprogress = conf.showprogress;

%% ===== Computation ====================================================
L = (Nse + 1).^2;
a = zeros(L,1);
b = zeros(L,1);
for n=0:Nse
  a_denum = sqrt((2*n+1)*(2*n+3));
  b_denum = sqrt((2*n-1)*(2*n+1));
  for m=0:n   
    l_plus = (n + 1).^2 - (n - m);
    l_minus = (n + 1).^2 - (n + m);    

    a(l_plus) = sqrt((n+1+m)*(n+1-m))./a_denum;
    a(l_minus)= a(l_plus);
    b(l_plus) = sqrt((n-m-1)*(n-m))./b_denum;
    b(l_minus) = -sqrt((n+m-1)*(n+m))./b_denum;
  end
  if showprogress, progress_bar(l_plus,L); end % progress bar
end