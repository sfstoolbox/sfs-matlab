function [ir,delay] = wfs_preequalization(ir,conf)
%WFS_PREEQUALIZATION applies a pre-equalization filter for WFS
%
%   Usage: ir = wfs_preequalization(ir,conf)
%
%   Input parameters:
%       ir      - signal to which the pre-equalization filter should be applied
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       ir      - signal with applied pre-equalization
%       delay   - additional delay added by pre-equalization / s
%
%
%   WFS_PREEQUALIZATION(ir,conf) applies the pre-equalization filter for
%   Wave Field Synthesis to the given impulse response.
%
%   See also: wfs_fir_prefilter, wfs_iir_prefilter, driving_function_imp_wfs

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


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargmatrix(ir);
isargstruct(conf);


%% ===== Configuration ==================================================
usehpre = conf.wfs.usehpre;
fs = conf.fs;


%% ===== Computation =====================================================
% Check if we should procide
if ~usehpre
    delay = 0;
    return;
end
% Store original length
len_ir = size(ir,1);
% Get the filter
if strcmp('FIR',conf.wfs.hpretype)
    % Get FIR filter
    hpre = wfs_fir_prefilter(conf);
    % Apply filter
    ir = convolution(hpre,ir);
    % Delay in s added by filter
    delay = conf.wfs.hpreFIRorder/2 / fs;
elseif strcmp('IIR',conf.wfs.hpretype)
    if len_ir == 1
        % Happens when called from driving_function_imp_wfs()
        % Zeropadding to length conf.N
       ir = [ir; zeros(conf.N-1,size(ir,2))];
    end
    % Get IIR filter
    hpre = wfs_iir_prefilter(conf);
    % Apply filter
    ir = sosfilt(hpre.sos,ir,1);
    % IIR is minimum phase, so no proper delay introduced
    delay = 0;
else
    error('%s: %s is an unknown filter type.',upper(mfilename),hpretype);
end
% Correct length of ir
if strcmp('FIR',conf.wfs.hpretype) && (len_ir>length(hpre)+1)
    ir = ir(1:len_ir,:);
end
