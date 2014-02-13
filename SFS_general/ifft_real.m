function y = ifft_real(x,N)
%IFFT_REAL computes an ifft for real value signals
%
%   Usage: y = ifft_real(x,samples)
%
%   Input parameters:
%       x       - matrix x with signals as columns
%       samples - desired length of the output signals
%
%   Output parameters:
%       y       - matrix y with signals as columns
%
%   IFFT_REAL(x,samples) computes the ifft of the real signals in x. The signals
%   have to be the columns of x.
%
%   see also: fft_real, convolution

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


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargmatrix(x);
isargpositivescalar(N);


%% ===== Computation =====================================================
% Force IFFT along dimension 1
if rem(N,2)==0
  y = [x; flipud(conj(x(2:end-1,:)))];
else
  y = [x; flipud(conj(x(2:end,:)))];
end;

y = real(ifft(y,N,1));


