function B = sphexp_mono_multiscatter(A, xq, R, sigma, f, conf)

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
nargmin = 5;
nargmax = 6;
narginchk(nargmin,nargmax);
isargmatrix(A,xq);
isargpositivescalar(f);
isargvector(R, sigma);
if nargin<nargmax
  conf = SFS_config;
else
  isargstruct(conf);
end
if length(R) == 1
  R = repmat(R,[1, size(A,2)]);
elseif length(R) ~= size(A,2)
  error('%s: Length of R does not match size of A',upper(mfilename));
end
if length(sigma) == 1
  sigma = repmat(sigma,[1, size(A,2)]);
elseif length(sigma) ~= size(A,2)
  error('%s: Length of R does not match size of A',upper(mfilename));
end

%% ===== Configuration ==================================================
Nse = conf.scattering.Nse;

%% ===== Computation ====================================================
if size(A,1) ~= (Nse + 1)^2
  error('%s: size of A does not match Nse',upper(mfilename));
end

Nnm = size(A,1);
Nq = size(A,2);
L = zeros(Nq*Nnm);

E = ones(Nnm,1);
for qdx=1:Nq
  selectq = ((qdx-1)*Nnm+1):(qdx*Nnm);
  L(selectq,selectq) = diag(1./sphexp_mono_scatter(E, R(qdx), sigma(qdx), f, conf));  
end

for qdx=1:Nq
  selectq = ((qdx-1)*Nnm+1):(qdx*Nnm);
  for pdx=(qdx+1):Nq
    selectp = ((pdx-1)*Nnm+1):(pdx*Nnm);    
    [SRpq, SRqp] = sphexp_mono_translation(xq(qdx,:)-xq(pdx,:), 'SR', f, conf);
    L(selectp,selectq) = -SRpq;
    L(selectq,selectp) = -SRqp;
  end
end

A = reshape(A,[],1);

B = L\A;
B = reshape(B,Nnm,Nq);
