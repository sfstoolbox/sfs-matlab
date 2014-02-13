function sig = fix_length(sig,N)
%FIX_LENGTH pads zeros or removes entries from the signal according to length N
%
%   Usage: sig = fix_length(sig,N)
%
%   Input parameters:
%       sig - input signal (matrix with sigs as columns)
%       N   - number of samples size(sig,1) should be
%
%   Output paramteres:
%       sig - corrected sig
%
%   FIX_LENGTH(sig,N) pads zeros or removes the end of the given signal in
%   order to have a sig with a size(sig,1)==N.
%
%   see also: fix_irs_length, convolution

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
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);


%% ===== Fix IR ==========================================================
% length of IR
samples = size(sig,1);
channels = size(sig,2);

if samples<N
    % append zeros if to short
    sig = [sig; zeros(N-samples,channels)];
else
    % remove the end of the IR, if to long
    sig = sig(1:N,:);
end
