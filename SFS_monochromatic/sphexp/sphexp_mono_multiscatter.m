function B = sphexp_mono_multiscatter(A, xq, R, sigma, f, conf)
%SPHEXP_MONO_MULTISCATTER compute the singular spherical expansion of a sphere-
%scattered sound field
%
%   Usage: B = sphexp_mono_multiscatter(A, xq, R, sigma, f, conf)
%
%   Input parameters:
%       A           - regular spherical expansion coefficients of incident 
%                     sound fields arriving at each sphere [n x Nq]
%       xq          - positions of the sphere / m [Nq x 3]
%       R           - radii of the spheres / m [Nq x 1]
%       sigma       - admittances of the sphere / m [Nq x 1]
%       f           - frequency / Hz
%
%   Output parameters:
%       B           - regular spherical expansion coefficients of sound fields
%                     scattered at each sphere [n x Nq]
%
%   SPHEXP_MONO_MULTISCATTER(A, xq, R, sigma, f, conf)

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
nargmin = 6;
nargmax = 6;
narginchk(nargmin,nargmax);
isargpositivescalar(f);
isargvector(R, sigma);
isargmatrix(A, xq);
isargstruct(conf);

L = size(A,1);
isargsquaredinteger(L);

if length(R) == 1
  R = repmat( R, [1, size(A,2)] );
elseif length(R) ~= size(A,2)
  error('%s: Length of R does not match size of A',upper(mfilename));
end

if length(sigma) == 1
  sigma = repmat( sigma, [1, size(A,2)] );
elseif length(sigma) ~= size(A,2)
  error('%s: Length of R does not match size of A',upper(mfilename));
end

%% ===== Computation ====================================================
Nq = size(A,2);
AB = zeros(Nq*L);

E = ones(L,1);
for qdx=1:Nq
  selectq = ((qdx-1)*L+1):(qdx*L);
  AB(selectq,selectq) = ...
    diag(1./sphexp_mono_scatter(E, R(qdx), sigma(qdx), f, conf));  
end

for qdx=1:Nq
  selectq = ((qdx-1)*L+1):(qdx*L);
  for pdx=(qdx+1):Nq
    selectp = ((pdx-1)*L+1):(pdx*L);    
    [SRpq, SRqp] = sphexp_mono_translation(xq(qdx,:)-xq(pdx,:), 'SR', f, conf);
    AB(selectp,selectq) = -SRpq;
    AB(selectq,selectp) = -SRqp;
  end
end

A = reshape(A,[],1);

B = AB\A;
B = reshape(B,L,Nq);
