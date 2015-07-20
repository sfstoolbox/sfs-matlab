function Pnm = sphexp_mono_nfchoa_sht(Dnm,mode,f,conf)
%SPHEXP_MONO_NFCHOA_SHT yields spherical expansion coefficients of a sound field
%resulting from of the nfchoa driving function given as a spherical harmonics
%transform
%
%   Usage: Pnm = sphexp_mono_nfchoa_sht(Dnm,mode,f,conf)
%
%   Input parameters:
%       Dnm         - spherical harmonics transform of nfchoa driving function
%       mode        - 'R' for regular, 'S' for singular
%       f           - frequency / Hz [Nf x 1] or [1 x Nf]
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       Pnm         - spherical expansion coefficients of a sound field
%                     reproduced by nfchoa driving function
%
%   SPHEXP_MONO_NFCHOA_SHT(Dnm,mode,f,conf)

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
isargmatrix(Dnm);
isargsquaredinteger(size(Dnm,1));
isargvector(f);
isargchar(mode);
if ~strcmp('R', mode) && ~strcmp('S', mode)
  error('%s: unknown mode (%s)!', upper(mfilename), mode);
end
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end

%% ===== Configuration ==================================================
r0 = conf.secondary_sources.size / 2;
dimension = conf.dimension;

%% ===== Variables ======================================================
Nse = sqrt(size(Dnm,1)) - 1;

%% ===== Computation ====================================================
if strcmp('2D',dimension)
  % === 2-Dimensional ==================================================

  error('%s: 2D not supported.',upper(mfilename));

elseif strcmp('2.5D',dimension)
  % === 2.5-Dimensional ================================================

  % regular/singular expansion of 3D Green's function
  Gnm = 2*pi*r0*sphexp_mono_ps([r0, 0, 0], mode, Nse, f, [0,0,0], conf);
elseif strcmp('3D',dimension)
  % === 3-Dimensional ==================================================

  % regular/singular expansion of 3D Green's function
  Gnm = sphexp_mono_ps([0, 0, r0], mode, Nse, f, [0,0,0], conf);

  for n=0:Nse
    v = sphexp_index(-n:n,n);
    w = sphexp_index(0,n);
    Gnm(v,:) = 2*pi*r0^2*sqrt(4*pi / (2*n+1))*Gnm(w,:);
  end
else
  error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end

Pnm = Gnm .* Dnm;

end