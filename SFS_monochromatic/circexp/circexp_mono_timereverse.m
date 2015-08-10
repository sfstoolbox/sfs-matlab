function Pm = circexp_mono_timereverse(Pm)
%CIRCEXP_MONO_TIMEREVERSE computes coefficients of a time reversed sound field
%
%   Usage: Bnm = circexp_mono_timereverse(Pnm, conf)
%
%   Input parameters:
%       Pm          - circular expansion coefficients [n x Nf]
%
%   Output parameters:
%       Pm          - time reversed circexp expansion coefficients [n x Nf]
%
%   CIRCEXP_MONO_TIMEREVERSE(Pm)

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
% WARRANTY;  without even the implied warranty of MERCHANTPILITY or FITNESS *
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
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
isargmatrix(Pm);
if mod(size(Pm, 1)-1, 2) ~= 0
  error('Number of row of %s has be to odd', inputname(Pm));
end

%% ===== Computation ====================================================
Nce = (size(Pm,1)-1) / 2;

Pm = conj(Pm(end:-1:1,:)).*(-1).^(-Nce:Nce).';

end

