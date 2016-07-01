function Pnm = sphexp_mono_cht(Dm, mode, f, conf)
%SPHEXP_MONO_CHT yields spherical expansion coefficients of a sound field
%resulting from of a driving function given as a circular harmonics transform
%
%   Usage: Pnm = sphexp_mono_cht(Dm, mode, f, conf)
%
%   Input parameters:
%       Dm          - circular harmonics transform of driving function
%       mode        - 'R' for regular, 'S' for singular
%       f           - frequency / Hz [Nf x 1] or [1 x Nf]
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       Pnm         - spherical expansion coefficients of a synthesized
%                     sound field reproduced by nfchoa driving function
%
%   SPHEXP_MONO_CHT(Dm, mode, f, conf) computes the spherical expansion 
%   coefficients of sound field reproduced by a circular secondary source
%   distribution consisting of SPHERICAL monopoles. The driving function is
%   given as its circular harmonics transform.

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
isargmatrix(Dm);
isargvector(f);
isargchar(mode);
if ~strcmp('R', mode) && ~strcmp('S', mode)
  error('%s: unknown mode (%s)!', upper(mfilename), mode);
end
isargstruct(conf);

%% ===== Configuration ==================================================
R0 = conf.secondary_sources.size / 2;

%% ===== Variables ======================================================
Nse = (size(Dm, 1)-1)/2;

%% ===== Computation ====================================================
% regular/singular spherical expansion of 3D Green's function
Gnm = sphexp_mono_ps([R0, 0, 0], mode, Nse, f, [0,0,0], conf);
% regular/singular spherical expansion of synthesised sound field
Pnm = zeros(size(Gnm));
for n=0:Nse
  v = sphexp_index(-n:n,n);
  Pnm(v,:) = 2*pi*R0*Gnm(v,:).*Dm((-n:n)+Nse+1,:);
end

end