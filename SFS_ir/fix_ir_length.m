function ir = fix_ir_length(ir,N)
%FIX_IR_LENGTH pads zeros or removes entries from the IR according to length N
%
%   Usage: ir = fix_ir_length(ir,N)
%
%   Input parameters:
%       ir  - impulse response (IR)
%       N   - number of samples the calculated BRIR should have
%
%   Output paramteres:
%       ir  - corrected IR
%
%   FIX_IR_LENGTH(IR,N) pads zeros or removes the end of the given IR in
%   order to have a IR with a length of N.
%
%   see also: fix_irs_length, get_ir, read_irs

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);


%% ===== Fix IR ==========================================================
% length of IR
samples = size(ir,1);
channels = size(ir,2);

if samples<N
    % append zeros if to short
    ir = [ir; zeros(N-samples,channels)];
else
    % remove the end of the IR, if to long
    ir = ir(1:N,:);
end
