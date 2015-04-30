function D = driving_function_mono_nfchoa_sphexp(x0, Pnm,f,conf)
%computes the nfchoa driving functions for a sound field expressed by regular 
%spherical expansion coefficients.
%
%   Usage: D = driving_function_mono_nfchoa_sphexp(x0,Pnm,f,N,conf)
%
%   Input parameters:
%       x0          - position of the secondary sources / m [nx3]
%       Pnm         - regular spherical expansion coefficients of sound field
%       f           - frequency in Hz
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_NFCHOA_SPHEXP(x0,Pnm,f,N,conf) returns NFCHOA
%   driving signals for the given secondary sources, the virtual sound expressed
%   by regular spherical expansion coefficients and the frequency f.
%
%   see also: driving_function_mono_wfs_sphexp

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
isargmatrix(x0);
isargvector(Pnm);
isargpositivescalar(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end
if mod(sqrt(size(Pnm, 1)),1) ~= 0
  error(['%s: number of columns of Pnm (%s) is not the square of an', ...
    'integer.'], upper(mfilename), sqrt(size(Pnm, 1)));
end

%% ===== Configuration ==================================================
c = conf.c;
dimension = conf.dimension;
Xc = conf.secondary_sources.center;

%% ===== Variables ======================================================
Nse = sqrt(size(Pnm, 1))-1;

%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

% secondary source positions
x00 = bsxfun(@minus,x0,Xc);
[phi0,theta0,r0] = cart2sph(x00(:,1),x00(:,2),x00(:,3));

% frequency depended stuff
omega = 2*pi*f;
k = omega/c;
kr0 = k.*r0;
% initialize empty driving signal
D = zeros(size(x0,1),1);

%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

if strcmp('2D',dimension)
  % === 2-Dimensional ==================================================
  to_be_implemented; 
  
elseif strcmp('2.5D',dimension)
  % === 2.5-Dimensional ================================================
  %                                m
  %                       __      P
  %               1      \         |m|
  % D(phi0,w) = ------   /__    -------- e^(i m (phi0))
  %             2pi r0  m=-N..N    m
  %                               G
  %                                |m|
  %
  % with regular spherical expansion of the sound field:
  %          \~~   N \~~   n   m           m
  % P(x,w) =  >       >       P  j (kr) . Y  (theta, phi)
  %          /__ n=0 /__ m=-n  n  n        n 
  %
  % and 3D free field Green's Function:
  %             \~~ oo  \~~   n   m           m
  % G  (x0,f) =  >       >       G  . j (kr) . Y  (theta, phi)
  %  ps         /__ n=0 /__ m=-n  n    n        n
  %
  % with the regular expansion coefficients (Gumerov2004, eq. 3.2.2):
  %    m               (2)       -m
  %   G  = -i  . k  . h   (kr0) Y  (pi/2, 0)
  %    n               n         n
  
  for m=-Nse:Nse
    D = D + sphexp_access(Pnm, m) ./ ...
      (sphbesselh(abs(m),2,kr0) .* sphharmonics(abs(m),-m, 0, 0) ) ...
      .* exp(1i.*m.*phi0);      
  end
  
  D = D./(-1j*2.*pi*kr0);  
  
elseif strcmp('3D',dimension)
  % === 3-Dimensional ==================================================
  
  to_be_implemented;  
else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
