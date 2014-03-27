function [SRpq, SRqp, T] = sphexpSR_mono(xpq, f, conf)
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
nargmin = 2;
nargmax = 3;
narginchk(nargmin,nargmax);
isargposition(xpq);
if nargin<nargmax
  conf = SFS_config;
else
  isargstruct(conf);
end

%% ===== Configuration ==================================================
showprogress = conf.showprogress;
Nse = conf.scattering.Nse;

%% ===== Computation ====================================================
% convert (xpq) into spherical coordinates
r = sqrt(sum(xpq.^2));
phi = atan2(xpq(2),xpq(1));
theta = asin(xpq(3)./r);

% frequency
k = 2*pi*f/conf.c;
kr = k*r;

L = (2*Nse + 1).^2;
S = zeros(L);
T = zeros(L);

[a, b] = sphexpSR_ab(2*Nse,conf);
%% ===== Sectorial Coefficients =========================================
% for n=0, m=0
for l=0:2*Nse
  Hn  = sqrt(4*pi)*sphbesselh(l,2,kr);
  for s=0:l
    Ynm = sphharmonics(l,-s, theta, phi);  % spherical harmonics
    % +s
    v = sphexp_index(s,l);
    S(v,1) = (-1).^s*Hn.*Ynm;
    T(v,1) = 4;
    % -s
    v = sphexp_index(-s,l);
    S(v,1) = (-1).^(-s)*Hn.*conj(Ynm);
    T(v,1) = 4;
  end
  if showprogress, progress_bar(v,L); end % progress bar
end

% for n=|m|
for m=0:Nse-1
  bm = -1./sphexp_access(b,-m-1);
  for l=0:(2*Nse-m-1)
    for s=-l:l
      % +m
      [v, w] = sphexp_index(s,l,m+1);
      S(v,w) = bm...
        *(sphexp_access(b,-s,l)     * sphexp_access(S,s-1,l-1,m)...
        + sphexp_access(b,s-1,l+1)  * sphexp_access(S,s-1,l+1,m));
      T(v,w) = 1;
      % -m
      [v, w] = sphexp_index(s,l,-m-1);
      S(v,w) = bm...
        *(sphexp_access(b,s,l)      * sphexp_access(S,s+1,l-1,-m)...
        + sphexp_access(b,-s-1,l+1) * sphexp_access(S,s+1,l+1,-m));
      T(v,w) = 1;
    end
  end
  if showprogress, progress_bar(m,Nse); end % progress bar
end
% symmetry relation: S(s, l, m, n) = (-1)^(l+n)*S(-m, n, -s, l)
for s=0:Nse
  l = s;
  for n=0:(2*Nse-s)
    for m=-n:n
      % +s
      [v, w] = sphexp_index( s, l, m, n);
      S(v,w) = (-1)^(l+n).*sphexp_access(S,-m,n,-s,l);
      T(v,w) = 2;
      % -s
      [v, w] = sphexp_index(-s, l, m, n);
      S(v,w) = (-1)^(l+n).*sphexp_access(S,-m,n, s,l);
      T(v,w) = 2;
    end
  end
end

%% ===== Tesseral Coefficients ==========================================
for m=-Nse:Nse
  for s=-Nse:Nse
    for n=abs(m):Nse-1
      amn1 = -1./sphexp_access(a,m,n);
      amn2 = sphexp_access(a,m,n-1);
      for l=(abs(s)+1):(2*Nse-n-1)
        % +m, +s
        [v, w] = sphexp_index(s,l, m, n+1);
        S(v,w) = amn1...
          * (-amn2                   * sphexp_access(S,s,l  ,m,n-1) ...
          +  sphexp_access(a,s,l)    * sphexp_access(S,s,l+1,m,n) ...
          -  sphexp_access(a,s,l-1)  * sphexp_access(S,s,l-1,m,n));
        T(v,w) = 3;       
      end
    end
  end  
  if showprogress, progress_bar(m+Nse,2*Nse); end % progress bar
end

%% ======================================================================
L = (Nse + 1).^2;
SRpq = S(1:L,1:L);  % SR(xq-xp)
SRqp = zeros(L);    % SR(xp-xq)
for n=0:Nse
  for l=0:Nse
    for m=-n:n
      for s=-l:l
        [v, w] = sphexp_index(s,l,m,n);
        SRqp(v,w) = (-1).^(l+n)*sphexp_access(SRpq,s,l,m,n);
      end
    end
  end
end


end


