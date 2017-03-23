function sig =  bandpass(sig,flow,fhigh,conf)
%BANDPASS filters a signal by a bandpass
%
%   Usage: sig = bandpass(sig,flow,fhigh,conf)
%
%   Input parameters:
%       sig    - input signal (matrix)
%       flow   - start frequency of bandpass
%       fhigh  - stop frequency of bandpass
%       conf   - configuration struct (see SFS_config)
%
%   Output parameters:
%       sig    - filtered signal
%
%   BANDPASS(sig,flow,fhigh,conf) filters the given signal with a bandpass
%   filter with cutoff frequencies of flow and fhigh.
%
%   See also: sound_field_imp_wfs

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargmatrix(sig);
isargpositivescalar(flow,fhigh);
isargstruct(conf);


%% ===== Configuration ==================================================
fs = conf.fs;
N = 128;


%% ===== Computation =====================================================
% Get frequency range
range = fhigh-flow;
% Calculate scaling factor for frequency range
scaling = range/20000;
% Design bandpass filter
Hf = [0 2*flow/fs (2+2*scaling)*flow/fs (2-0.2*scaling)*fhigh/fs 2*fhigh/fs 1];
Hm = [0 0 1 1 0 0];
b = fir2(N,Hf,Hm);
% Filter signal
sig = convolution(sig,b);
% Compensate for delay & truncate result
sig = sig(N/2:end-(N/2)-1,:);
