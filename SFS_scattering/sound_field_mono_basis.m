function P = sound_field_mono_basis(AB, JH2n, Y,conf)
%SOUND_FIELD_MONO_BASIS simulates a sound field with cylindrical/spherical
%basic functions
%
%   Usage: P = sound_field_mono_basis(AB, JH2n, Y,conf)
%
%   Input parameters:
%       AB          - regular/singular cylindrical/spherical expansion coefficients
%       JH2n        - cell array of cylindrical/spherical bessel/hankel(2nd kind) functions
%       Y           - cell array of cylindrical/spherical harmonics
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       P           - resulting soundfield
%
%   SOUND_FIELD_MONO_BASIS(AB, JH2n, Y,conf)
%
%   see also: cylbasis_mono_XYZgrid, sphbasis_mono_XYZgrid

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
isargvector(AB);
if nargin<nargmax
  conf = SFS_config;
else
  isargstruct(conf);
end
if length(AB) ~= length(Y)
  error('%s: length missmatch.',upper(mfilename));
end

%% ===== Configuration ==================================================
% Plotting result
showprogress = conf.showprogress;

%% ===== Computation ====================================================
L = length(Y);
N = length(JH2n);
P = zeros(size(JH2n{1}));

if L == N
  % cylindrical basic functions
  for l=1:L    
    P = P + AB(l)*(JH2n{l}.*Y{l});
    if showprogress, progress_bar(l,L); end  % progress bar
  end  
elseif L == N^2
  % spherical basic functions
  l = 0;
  for n=0:N-1
    for m=-n:n
      l=l+1;
      if showprogress, progress_bar(l,L); end  % progress bar
      P = P + AB(l)*(JH2n{n+1}.*Y{l});
    end
  end  
else
  error(['%s: lengths of arrays are not suited for cylindrical or .' ...
    , 'spherical basic funcions'],upper(mfilename));
end

end

