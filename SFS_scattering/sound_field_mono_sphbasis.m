function P = sound_field_mono_sphbasis(ABnm, jh2n, Ynm,conf)
%SOUND_FIELD_MONO_BASIS simulates a sound field with cylindrical/spherical
%basic functions
%
%   Usage: P = sound_field_mono_sphbasis(AB, jh2n, Y,conf)
%
%   Input parameters:
%       ABnm        - regular/singular spherical expansion coefficients
%       jh2n        - cell array of spherical bessel/hankel(2nd kind) functions
%       Ynm         - cell array of spherical harmonics
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       P           - resulting soundfield
%
%   SOUND_FIELD_MONO_SPHBASIS(AB, JH2n, Y,conf)
%
%   see also: sphbasis_mono_grid

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
isargvector(ABnm);
if nargin<nargmax
  conf = SFS_config;
else
  isargstruct(conf);
end
isargequallength(ABnm, Ynm);
L = length(Ynm);
Nse = length(jh2n) - 1;
if (L ~= (Nse + 1)^2)
  error('%s, length(Y) has to be (length(jhsn)+1).^2!',upper(mfilename));
end

%% ===== Configuration ==================================================
% Plotting result
showprogress = conf.showprogress;
usenormalisation = conf.usenormalisation;

%% ===== Computation ====================================================
P = zeros(size(jh2n{1}));

% spherical basic functions
l = 0;
for n=1:Nse
  for m=-n:n
    l=l+1;
    if showprogress, progress_bar(l,L); end  % progress bar
    P = P + ABnm(l)*(jh2n{n}.*Ynm{l});
  end
end

if usenormalisation
  P = P./max(abs(P(:)));
end

end

