function ir = fix_ir_length(ir,N,dt)
%FIX_IR_LENGTH pads zeros or removes entries from the IR according to length N
%
%   Usage: ir = fix_ir_length(ir,N,[dt])
%
%   Input parameters:
%       ir  - impulse response (IR)
%       N   - number of samples the calculated BRIR should have
%       dt  - time delay for the given setup the IR will be shifted with, 
%             default: 0
%
%   Output paramteres:
%       ir  - corrected IR
%
%   FIX_IR_LENGTH(IR,N,DT) pads zeros or removes the end of the given IR in
%   order to have a IR with the correct length to feed it with the time delay dt
%   in the desired BRIR with the length N.
%
%   see also: brs_point_source, get_ir

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
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
nargmax = 3;
narginchk(nargmin,nargmax);
if nargin==nargmax-1
    dt = 0;
end


%% ===== Fix IR ==========================================================
% length of IR
lenir = length(ir(:,1));
% ensure integer delays
dt = round(dt);
% append zeros if to short
if(lenir<N+abs(dt))
    ir = cat(1,ir,zeros(N-lenir,2));
% remove the end of the IR, if to long
else
    ir=ir(1:N,:);
end
