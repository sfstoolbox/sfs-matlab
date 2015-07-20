function [a, b] = sphexp_translation_auxiliary(Nse,conf)
%Auxiliary coefficients for the recursive calculation of spherical translation coefficients
%
%   Usage: [a, b] = sphexp_translation_auxiliary(Nse,conf)
%
%   Input parameters:
%       Nse         - maximum degree of the auxiliary functions (optional)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       a           - auxiliary coefficients for the calculation of tesseral
%                     spherical translation coefficients
%       b           - auxiliary coefficients for the calculation of sectorial
%                     spherical translation coefficients
%
%   SPHEXP_TRANSLATION_AUXILIARY(Nse,conf)
%
%   Auxiliary coefficients for the calculation of tesseral spherical
%   translation coefficients (Gumerov2004, eq. 2.2.8 and 3.2.65):
%
%         +------------------+
%    m    |(n+1+|m|)(n+1-|m|)           _
%   a  =  |------------------  for |m|  < n
%    n   \|   (2n+1)(2n+3)
%
%         0                    else
%
%   Auxiliary coefficients for the calculation of sectorial spherical
%   translation coefficients (Gumerov2004, eq. 2.2.10 and 3.2.78/79):
%
%          +------------+
%          |(n-m)(n-m-1)              _    _
%          |------------       for 0  < m  < n
%         \|(2n-1)(2n+1)
%
%          +------------+
%    m     |(n-m)(n-m-1)               _
%   b  = - |------------       for -n  < m  < 0
%    n    \|(2n-1)(2n+1)
%
%          0                   else
%
%   The coefficients are stored in linear arrays with index l resulting from
%   m and n:
%
%         m         m               2
%   a  = a  ; b  = a  with l = (n+1)  - (n - m)
%    l    n    l    n
%
%   References:
%       Gumerov,Duraiswami (2004) - "Fast Multipole Methods for the
%                                    Helmholtz Equation in three
%                                    Dimensions", ELSEVIER
%
%   see also: sphexp_mono_ps, sphexp_mono_pw
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
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
isargpositivescalar(Nse);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end

%% ===== Configuration ==================================================
showprogress = conf.showprogress;

%% ===== Computation ====================================================
L = (Nse + 1).^2;
a = zeros(L,1);
b = zeros(L,1);
for n=0:Nse
  a_denum = (2*n+1)*(2*n+3);
  b_denum = (2*n-1)*(2*n+1);

  m = -n:n;
  v = sphexp_index(m,n);

  a(v) = sqrt( (n+1+abs(m)).*(n+1-abs(m))./a_denum );
  b(v) = ( 1-(m<0)*2 ) .* sqrt( (n-m-1).*(n-m)./b_denum );

  if showprogress, progress_bar(v(end),L); end % progress bar
end
