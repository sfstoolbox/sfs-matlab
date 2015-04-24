function [EF, EFm] = sphexp_mono_translation(t, mode, f, conf)
% Spherical translation coefficients (multipole re-expansion)
%
%   Usage: [EF, EFm] = sphexp_mono_translation(t, mode, f, conf)
%
%   Input parameters:
%       t           - translatory shift [1x3] / m                    
%       mode        - 'RS' for regular-to-singular reexpansion
%                     'RR' for regular-to-regular reexpansion
%                     'SR' for singular-to-regular reexpansion
%                     'SS' for singular-to-singular reexpansion
%       f           - frequency / Hz
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       EF          - singular spherical expansion coefficients of
%                     scattered field
%  
%  SPHEXP_MONO_TRANSLATION(t, mode, f, conf) computes the spherical re-expansion
%  coefficients to perform as translatory shift of spherical basis function.
%  Multipole Re-expansion computes the spherical basis function for a shifted
%  coordinate system (x+t) based on the original basis functions for (x). 
%
%   m          \~~ inf \~~ l         s,m     s
%  E (x + t) =  >       >       (E|F)   (t) F (x)
%   n          /__ l=0 /__ s=-l      l,n     l
%
%  where {E,F} = {R,S}. R denotes the regular spherical basis function, while
%  S symbolizes the singular spherical basis function. Note that (S|S) and 
%  (S|R) are respectively equivalent to (R|R) and (R|S). 
%  The reexpansion coefficients can seperated into sectorial 
%  (n = abs|m| and/or l = abs|s|) and tesseral (else) coefficients. Latter will
%  only be computed, if conf.dimensions == '3D'.
%
%  References:
%     Gumerov,Duraiswami (2004) - "Fast Multipole Methods for the 
%                                   Helmholtz Equation in three 
%                                   Dimensions", ELSEVIER
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
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargposition(t);
isargchar(mode);
isargpositivescalar(f);
if nargin<nargmax
  conf = SFS_config;
else
  isargstruct(conf);
end

%% ===== Configuration ==================================================
showprogress = conf.showprogress;
Nse = conf.scattering.Nse;
dimension = conf.dimension;

%% ===== Variables ======================================================
% convert (xpq) into spherical coordinates
r = sqrt(sum(t.^2));
phi = atan2(t(2),t(1));
theta = asin(t(3)./r);

% frequency
k = 2*pi*f/conf.c;
kr = k*r;

L = (2*Nse + 1).^2;
S = zeros(L);

% Auxilary Coefficients for Computations
[a, b] = sphexp_translation_auxiliary(2*Nse,conf);

% select suitable basis function
if strcmp('RR', mode) || strcmp('SS', mode)
  sphbasis = @sphbesselj;
elseif strcmp('SR', mode) || strcmp('RS', mode)
  sphbasis = @(nu,z) sphbesselh(nu,2,z);
else
  error('unknown mode:');
end

%% ===== Sectorial Coefficients =========================================

% for n=0, m=0 (Gumerov2004, eq. 3.2.5)
%
%      s, 0        s, 0     ___     -l  -s
% (S|R)     = (R|S)     = \|4pi (-1)   S  (t)
%      l, 0        l, 0                 l
%
%      s, 0        s, 0     ___     -l  -s
% (S|S)     = (R|R)     = \|4pi (-1)   R  (t)
%      l, 0        l, 0                 l
%
% for s=0, l=0
%
%      0, n        0, m     ___  m
% (S|R)     = (R|S)     = \|4pi S  (t)
%      0, m        0, n          n
%
%      0, n        0, m     ___   m
% (S|S)     = (R|R)     = \|4pi R  (t)
%      0, m        0, n          n
%
for l=0:2*Nse  % n=l
  Hn  = sqrt(4*pi)*sphbasis(l,kr);  % radial basis function (see Variables)
  for s=0:l  % m=s
    % spherical harmonics: conj(Y_n^m) = Y_n^-m (Gumerov2004, eq. 2.1.59)
    Ynm = sphharmonics(l,s, theta, phi);

    v = sphexp_index(s,l);
    S(v,1) = (-1).^l*Hn.*Ynm;  % s,l, m=0, n=0
    S(1,v) = Hn.*conj(Ynm);   % s=0,l=0, m, n

    v = sphexp_index(-s,l);
    S(v,1) = (-1).^l*Hn.*conj(Ynm); % -s,l, m=0, n=0
    S(1,v) = Hn.*Ynm;  % s=0,l=0, -m, n  
  end
  if showprogress, progress_bar(v,L); end % progress bar
end

% for n=|m| (Gumerov2004, eq. 3.2.78)
%
%  -m-1      s,m+1      -s      s-1,m      s-1      s-1,m
% b     (E|F)(t)     = b   (E|F)(t)     - b    (E|F)(t)
%  m+1       l,|m+1|    l       l-1,|m|    l+1      l+1,|m|
%
%  -m-1      s,-m-1     s       s+1,-m     -s-1     s+1,-m
% b     (E|F)(t)     = b   (E|F)(t)     - b    (E|F)(t)
%  m+1       l,|m+1|    l       l-1,|m|    l+1      l+1,|m|
% 
% while {E,F} = {R,S}
%
for m=0:Nse-1
  % NOTE: sphexp_access(b, -m-1) == sphexp_access(b, -m-1, m+1)
  bm = 1./sphexp_access(b,-m-1); 
  for l=1:(2*Nse-m-1)
    for s=-l:l
      % +m
      [v, w] = sphexp_index(s,l,m+1);
      S(v,w) = bm * ( ...
        sphexp_access(b,-s ,l)    * sphexp_access(S,s-1,l-1,m) - ...
        sphexp_access(b,s-1,l+1)  * sphexp_access(S,s-1,l+1,m)...
        );  
      % -m
      [v, w] = sphexp_index(s,l,-m-1);
      S(v,w) = bm * ( ...
        sphexp_access(b,s   ,l)   * sphexp_access(S,s+1,l-1,-m) - ...
        sphexp_access(b,-s-1,l+1) * sphexp_access(S,s+1,l+1,-m)...
        );     
    end
  end
  if showprogress, progress_bar(m,Nse); end % progress bar
end

% for l=|s| using symmetry relation (Gumerov2004, eq. 3.2.49
%
%      s,m         |s|+n     -m,-s
% (E|F)(t)   = (-1)     (E|F)(t)
%      |s|,n                 n,|s|
%
% while {E,F} = {R,S}
%
for l=1:Nse
  for s=[-l,l]
    for n=1:(2*Nse-l)
      for m=(-n+1):(n-1)
        % +s
        [v, w] = sphexp_index( s, l, m, n);
        S(v,w) = (-1).^(l+n)*sphexp_access(S,-m,n,-s);
      end
    end
  end
end

%% ===== Tesseral Coefficients ==========================================
if strcmp('3D',dimension)

  for m=-Nse:Nse
    for s=-Nse:Nse
      
      lowerbound = Nse - max(abs(m),abs(s));
      % left propagation
      for n=(0:lowerbound-1)+abs(m)
        amn1 = 1./sphexp_access(a,m,n);
        amn2 = sphexp_access(a,m,n-1);        
        for l=(abs(s)-abs(m)+n+1):(2*Nse-n-1)
          [v, w] = sphexp_index(m, n+1, s, l);
          S(v,w) = amn1 * ( ...
            amn2                    * sphexp_access(S,m,n-1,s,l) - ...
            sphexp_access(a,s,l)    * sphexp_access(S,m,n  ,s,l+1)   + ...
            sphexp_access(a,s,l-1)  * sphexp_access(S,m,n  ,s,l-1) ...
          );
        end
      end
    end
    % up propagation
    for l=(0:lowerbound-1)+abs(s)
      asl1 = 1./sphexp_access(a,s,l);
      asl2 = sphexp_access(a,s,l-1);        
      for n=(abs(m)-abs(s)+l+2):(2*Nse-l-1)
        [v, w] = sphexp_index(s,l+1, m, n);
        S(v,w) = asl1 * ( ...
          asl2                    * sphexp_access(S,s,l-1,m,n)   - ...
          sphexp_access(a,m,n)    * sphexp_access(S,s,l  ,m,n+1) + ...
          sphexp_access(a,m,n-1)  * sphexp_access(S,s,l  ,m,n-1)...
          );
      end
    end  
    if showprogress, progress_bar(m+Nse,2*Nse); end % progress bar
  end

end
%% ====== Final Calculation Steps =======================================
L = (Nse + 1)^2;
EF = S(1:L,1:L);  % (E|F)(t)
% SR(-t)
EFm = zeros(L);  
for n=0:Nse
  for l=0:Nse
    for m=-n:n
      for s=-l:l
        [v, w] = sphexp_index(s,l,m,n);
        EFm(v,w) = (-1).^(l+n)*sphexp_access(EF,s,l,m,n);
      end
    end
  end
end

end
