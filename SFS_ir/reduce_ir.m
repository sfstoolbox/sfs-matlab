function short_ir = reduce_ir(ir,fs,nsamples,conf)
%REDUCE_IR resamples and shortens a IR
%
%   Usage: ir = reduce_ir(ir,fs,nsamples,conf)
%
%   Input parameters:
%       ir          - two channel impulse response signal
%       fs          - sampling rate of the target impulse response / Hz
%       nsamples    - length of the target impulse response
%       conf        - configuration struct (see SFS_config)
%
%   Output paramteres:
%       ir          - two channel impulse response signal
%
%   REDUCE_IR(ir,fs,nsamples,conf) shortens and resamples a given impulse
%   response. This can be useful for mobile phones.
%
%   See also: get_ir, shorten_ir

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
isargpositivescalar(fs,nsamples);
if ~isnumeric(ir) || size(ir,2)~=2
    error('%s: ir has to be an impulse response with samples x 2 size.', ...
        upper(mfilename));
end
isargstruct(conf);


%% ===== Configuration ==================================================
ofs = conf.fs;  % original fs
useplot = conf.plot.useplot;


%% ===== Computation ====================================================
% Resample HRIR
if ofs~=fs
    resamp_ir(:,1) = resample(ir(:,1),fs,ofs);
    resamp_ir(:,2) = resample(ir(:,2),fs,ofs);
else
    resamp_ir = ir;
end
% Window HRIR
win = hann_window(ceil(0.15*nsamples),ceil(0.10*nsamples),nsamples).^2;
% Find maximum of resampled HRIR
% Find maximum in each channel and calculate the mean of the index
[~,idx1] = max(abs(resamp_ir(:,1)));
[~,idx2] = max(abs(resamp_ir(:,2)));
idx = round((idx1+idx2)/2);
% Cut the HRIR around the maximum
% Leading zeros before idx
offset = 24;
short_ir(1:nsamples,1) = ...
    resamp_ir(idx-offset:idx+nsamples-offset-1,1) .* win;
short_ir(1:nsamples,2) = ...
    resamp_ir(idx-offset:idx+nsamples-offset-1,2) .* win;


%% ===== Plotting =======================================================
if(useplot)
    figure
    plot(resamp_ir(:,1),'-b'); hold on;
    plot(resamp_ir(:,2),'-r');
    figure
    plot(short_ir(:,1),'-b'); hold on;
    plot(short_ir(:,2),'r-');
end
