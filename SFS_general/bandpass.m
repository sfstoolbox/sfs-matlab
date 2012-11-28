function sig =  bandpass(sig,conf)
%BANDPASS filters a signal by a bandpass
%
%   Usage: sig = bandpass(sig,[conf])
%
%   Input parameters:
%       sig  - input signal (vector)
%       conf - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       sig  - filtered signal
%
%   BANDPASS(sig) filters the given signal with a bandpass filter
%   with cutoff frequencies of 10Hz and 20kHz.
%
%   see also: wave_field_imp_wfs_25d

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
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
isargvector(sig);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
fs = conf.fs;
N = 128;


%% ===== Computation =====================================================
% design bandpass filter
Hf = [0 2*10/fs 2*20/fs 2*18000/fs 2*20000/fs 1];
Hm = [0 0 1 1 0 0];
b = fir2(N,Hf,Hm);
% filter signal
sig = conv(sig,b);
% compensate for delay & truncate result
sig = sig(N/2:end-(N/2)-1);
