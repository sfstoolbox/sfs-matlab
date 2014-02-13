function z = convolution(x,y)
%CONVOLUTION convolve the signals x and y
%
%   Usage: z = convolution(x,y)
%
%   Input parameters:
%       x       - matrix/vector with signals as columns
%       y       - matrix/vector with signals as columns, note that only one of
%                 the signals can be a matrix
%
%   Output parameters:
%       z       - convolved signal
%
%   CONVOLUTION(x,y) convolves the signals given with x and y. One of the input
%   signals can be a matrix containing the signals as column vectors, the other
%   one has to be a column vector. The convolution is done in the frequency
%   domain and it is checked if we have only real signals to speed up the
%   calculation. The length of z is length(x)+length(y)-1.
%
%   see also: fft_real, ifft_real, fft, ifft

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


%% ===== Checking input parameters =======================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargmatrix(x,y);
% check if only one of the inputs is a matrix
if all(size(x)>1) && all(size(y)>1)
    error('%s: Only one of the inputs can be multi-dimensional.', ...
        upper(mfilename));
end
% ensure column vectors
if ~all(size(x)>1), x=column_vector(x); end
if ~all(size(y)>1), y=column_vector(y); end


%% ===== Computation =====================================================
% if one of the input signals is a matrix repmat the vector of the other signal
if all(size(x)>1)
    y = repmat(y,1,size(x,2));
elseif all(size(y)>1)
    x = repmat(x,1,size(y,2));
end
% length of output signal
N = size(x,1)+size(y,1)-1;
% convolve the signals in frequency domain
if isreal(x) && isreal(y)
    z = ifft_real(fft_real(fix_length(x,N)).*fft_real(fix_length(y,N)),N);
else
    z = ifft(fft(fix_length(x,N)).*fft(fix_length(y,N)));
end
