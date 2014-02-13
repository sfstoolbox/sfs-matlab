function sig =  bandpass(sig,flow,fhigh,conf)
%BANDPASS filters a signal by a bandpass
%
%   Usage: sig = bandpass(sig,flow,fhigh,[conf])
%
%   Input parameters:
%       sig    - input signal (matrix)
%       flow   - start frequency of bandpass
%       fhigh  - stop frequency of bandpass
%       conf   - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       sig    - filtered signal
%
%   BANDPASS(sig,flow,fhigh) filters the given signal with a bandpass filter
%   with cutoff frequencies of flow and fhigh.
%
%   see also: sound_field_imp_wfs

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
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargmatrix(sig);
isargpositivescalar(flow,fhigh);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
fs = conf.fs;
N = 128;


%% ===== Computation =====================================================
% get frequency range
range = fhigh-flow;
% calculate scaling factor for frequency range
scaling = range/20000;
% design bandpass filter
Hf = [0 2*flow/fs (2+2*scaling)*flow/fs (2-0.2*scaling)*fhigh/fs 2*fhigh/fs 1];
Hm = [0 0 1 1 0 0];
b = fir2(N,Hf,Hm);
% filter signal
sig = convolution(sig,b);
% compensate for delay & truncate result
sig = sig(N/2:end-(N/2)-1,:);
