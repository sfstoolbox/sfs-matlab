function P = sound_field_mono_circbasis(Pm, JH2m, Ym)
%simulates a sound field with spherical basis functions
%
%   Usage: P = sound_field_mono_circbasis(Pm, JH2m, Ym)
%
%   Input parameters:
%       Pm          - regular/singular circular expansion coefficients
%       JH2m        - cell array of cylindrical bessel/hankel(2nd kind) functions
%       Ym          - cell array of cylindrical harmonics (exponential
%                     functions)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       P           - resulting soundfield
%
%   SOUND_FIELD_MONO_CIRCBASIS(Pm, JH2m, Ym)
%
%   see also: circbasis_mono_grid sound_field_mono_circexp

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
nargmax = 3;
narginchk(nargmin,nargmax);
isargvector(Pm);
isargequallength(Pm, Ym);
if length(Ym) ~= length(JH2m)
  error('%s, length(%s) has to be length(%s)!', upper(mfilename), ...
    inputname(2), inputname(3));
end

%% ===== Variables ======================================================
Nce = length(JH2m) - 1;

%% ===== Computation ====================================================
P = zeros(size(JH2m{1}));

l = 0;
for n=0:Nce
  l=l+1;    
  P = P + Pm(l).*(JH2m{l}.*Ym{l});
end

end