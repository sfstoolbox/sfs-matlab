function Dnm = driving_function_mono_nfchoa_sht_sphexp(Pnm, f, conf)
%DRIVING_FUNCTION_MONO_NFCHOA_SHT_SPHEXP computes the spherical harmonics 
%transform of nfchoa driving functions for a sound field expressed by regular 
%spherical expansion coefficients.
%
%   Usage: D = driving_function_mono_nfchoa_sht_sphexp(Pnm, f, conf)
%
%   Input parameters:
%       Pnm         - regular spherical expansion coefficients of virtual
%                     sound field [nxm]
%       f           - frequency in Hz [m x 1] or [1 x m]
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       Dnm         - regular spherical harmonics transform of driving
%                     function signal [n x m]
%
%   DRIVING_FUNCTION_MONO_NFCHOA_SHT_SPHEXP(Pnm, f, conf) returns spherical
%   harmonics transform of the NFCHOA driving function for a virtual sound 
%   expressed by regular spherical expansion coefficients and the frequency f.

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
nargmin = 3;
nargmax = 3;
narginchk(nargmin,nargmax);
isargmatrix(Pnm);
isargsquaredinteger(size(Pnm,1));
isargvector(f);
isargstruct(conf);
if size(Pnm,2) ~= length(f)
  error( '%s:number of rows in %s have to match length of %s', ...
    upper(mfilename), inputname(1), inputname(2) );
end

%% ===== Configuration ==================================================
c = conf.c;
dimension = conf.dimension;
driving_functions = conf.driving_functions;

Xc = conf.secondary_sources.center;
r0 = conf.secondary_sources.size / 2;
xref = conf.xref - Xc;
rref = norm( xref );  % reference radius
thetaref = asin( xref(3)./rref);  % reference elevation angle

%% ===== Variables ======================================================
Nse = sqrt(size(Pnm,1))-1;

%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

% frequency depended stuff
omega = 2*pi*row_vector(f);
k = omega./c;
kr0 = k.*r0;
krref = k.*rref;

%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

% initialize empty driving signal
Dnm = zeros(size(Pnm));

if strcmp('2.5D',dimension)
  % === 2.5-Dimensional ================================================
  
  if strcmp('default',driving_functions)
    
    if (rref == 0)
      % --- Xref == Xc -------------------------------------------------
      %                 m
      %                P
      %  n     1        |m|
      % D  = ------  --------
      %  m   2pi r0     m
      %                G
      %                 |m|
      %
      % with the regular expansion coefficients (Gumerov2004, eq. 3.2.2):
      %    m               (2)       -m
      %   G  = -i  . k  . h   (kr0) Y  (pi/2, 0)
      %    n               n         n
      
      for m=-Nse:Nse      
        v = sphexp_index(m, abs(m):Nse);  % v(1) contains index for n=abs(m)
      
        Dnm(v,:) = repmat( Pnm(v(1),:) ./ ( -1i.*k.* ...
          sphbesselh(abs(m),2,kr0).*sphharmonics(abs(m),-m, 0, 0) ) ...
          , length(v), 1);
      end      
    else
      % --- Xref ~= Xc --------------------------------------------------
      %
      % if the reference position is not in the middle of the array, 
      % things get a 'little' more complicated
      
      hn = zeros(size(Pnm));
      jn = hn;
      for n=0:Nse
        hn(n+1,:) = sphbesselh(n,2,kr0);
        jn(n+1,:) = sphbesselj(n,krref);
      end
    
      for m=-Nse:Nse
        Pm = zeros(1,length(f));
        Gm = Pm;
        % for theta=0 the legendre polynom is zero if n+m is odd
        for n=abs(m):2:Nse
          factor = jn(n+1,:) .* ...
            (-1).^(m) .* ...
            sqrt( (2*n+1) ./ (4*pi) ) .* ...
            sqrt( factorial(n-abs(m)) ./ factorial(n+abs(m)) ) .* ...
            asslegendre(n,abs(m), sin(thetaref));

          v = sphexp_index(m,n);
        
          Pm = Pm + Pnm(v,:) .* factor;

          Gm = Gm + (-1i*k) .* hn(n+1,:) .* sphharmonics(n, -m, 0, 0) ...
              .* factor;
        end

        v = sphexp_index(m, abs(m):Nse);        
        Dnm(v,:) = repmat( Pm./Gm, length(v), 1 );        
      end     
    end
    Dnm = Dnm./(2*pi*r0);
    
  else
    error('%s: %s, this type of driving function is not implemented ', ...
      upper(mfilename), driving_functions);
  end
  
elseif strcmp('3D',dimension)
  % === 3-Dimensional ==================================================

  if strcmp('default',driving_functions)
    
    for n=0:Nse
      Gn0 = (-1i*k) .* sphbesselh(n,2,kr0) .* sphharmonics(n, 0, pi/2, 0); 
      
      v = sphexp_index(-n:n, n);
      Dnm(v,:) = sqrt( (2*n+1) ./ (4*pi) ) .* Pnm(v,:) ./ ...
        repmat(Gn0,length(v),1);
    end    
    
    Dnm = Dnm./(2*pi*r0.^2);
  else
    error('%s: %s, this type of driving function is not implemented ', ...
      upper(mfilename), driving_functions);
  end
else
  error('%s: the dimension %s is unknown or not supported.',  ...
    upper(mfilename),dimension);
end
