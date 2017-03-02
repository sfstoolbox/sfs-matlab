function Dnm = driving_function_mono_wfs_sht_sphexp(Pnm, mode, f, conf)
%DRIVING_FUNCTION_MONO_WFS_SHT_SPHEXP computes spherical harmonics transform 
%the wfs driving functions for a sound field expressed by spherical expansion 
%coefficients.
%
%   Usage: Dnm = driving_function_mono_wfs_sht_sphexp(Pnm, mode, f, conf)
%
%   Input parameters:
%       Pnm         - regular/singular spherical expansion coefficients of 
%                     sound field
%       mode        - 'R' for regular expansion, 'S' for singular expansion
%       f           - frequency / Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       Dnm         - regular spherical harmonics transform of driving
%                     function signal [n x m]
%
%   DRIVING_FUNCTION_MONO_WFS_SHT_SPHEXP(Pnm, mode, f, conf)

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
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargmatrix(Pnm);
isargvector(f);
isargchar(mode);
isargstruct(conf);

%% ===== Configuration ==================================================
c = conf.c;
driving_functions = conf.driving_functions;
r0 = conf.secondary_sources.size / 2;

%% ===== Variables ======================================================
Nse = sqrt(size(Pnm, 1))-1;

% frequency depended stuff
k = 2.*pi.*f(:)./c;
kr0 = k.*r0;

% select suitable basis function
if strcmp('R', mode)
  sphbasis_derived = @(nu) k.*sphbesselj_derived(nu, kr0);
elseif strcmp('S', mode)
  sphbasis_derived = @(nu) k.*sphbesselh_derived(nu,2,kr0);
else
  error('%s: unknown mode (%s)!', upper(mfilename), mode);
end

%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

if (strcmp('default',driving_functions))
  % --- SFS Toolbox ------------------------------------------------
  %
  %             d    
  % D(x0, w) = ---- Ps(x0,w)
  %            d r0
  % with regular/singular spherical expansion of the sound field:
  %          \~~   N \~~   n   m  m
  % P(x,w) =  >       >       B  F (x) 
  %          /__ n=0 /__ m=-n  n  n
  %
  % where F = {R,S}.
  %
  % regular spherical basis functions:
  %  m                  m
  % R  (x) = j (kr)  . Y  (theta, phi)
  %  n        n         n     
  % singular spherical basis functions:
  %  m        (2)         m
  % S  (x) = h   (kr)  . Y  (theta, phi)
  %  n        n           n    

  Dnm = zeros( size(Pnm) );  
  for n=0:Nse
    v = sphexp_index(-n:n, n);  % m=-n:n
    
    Dnm(v,:) = Pnm(v,:) .* repmat( sphbasis_derived(n), 2*n+1, 1);
  end
  Dnm = -2.*Dnm;
  
else
  error(['%s: %s, this type of driving function is not implemented ', ...
    'for a 2.5D point source.'], upper(mfilename), driving_functions);
end
