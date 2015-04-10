function h = general_least_squares(samples,fractional_delay,passband_edge)
%GENERAL_LEAST_SQUARES fractional delay approximation using the general least
%squares method
%
%   Usage: h = general_least_squares(samples,fractional_delay,passband_edge)
%
%   Input parameters:
%       samples          - filter length (filter order N = samples-1)
%       fractional_delay - fractional delay (0 < x <= 1)
%       passband_edge    - passband edge of approximation (in [0 1])
%
%   Output parameters:
%       h                - filter coefficient vector h(1)...h(samples)
%
%   GENERAL_LEAST_SQUARES(samples,fractional_delay,passband_edge) calculates the
%   filter coefficients needed for the least squares fractional delay method
%   which is used in delayline() if conf.usefracdelay == true and
%   conf.fracdelay_method == 'least_squares'.
%
%   See also: delayline

% This code is based on hgls2() from Timo Laakso

%*****************************************************************************
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
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
% Disabled for time constrains
%nargmin = 3;
%nargmax = 3;
%narginchk(nargmin,nargmax);


%% ===== Computation =====================================================
N = samples-1;                   % filter order
M = N/2;                         % middle value
if (M-round(M))==0
    D = fractional_delay + M;    % integer part closest to middle
else
    D = fractional_delay + M-0.5;
end

cT = zeros(N+1,1);
p1 = cT;
cT(1)=passband_edge;
if round(D)==D
    p1(1) = passband_edge;
else
  p1(1) = ( sin(D*passband_edge*pi) )/(D*pi);
end
for k=1:N           % compute the elements of the Toeplitz matrix (vector)
  k1 = k+1;
  kD = k-D;
  cT(k1) = ( sin(k*passband_edge*pi) )/(k*pi);
  p1(k1) = ( sin(kD*passband_edge*pi) )/(kD*pi);
end
P = toeplitz(cT);
h = P\p1;
