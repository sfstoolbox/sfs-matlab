function b = pm_filter(order,wpass,wstop)
%PM_FILTER computes an FIR lowpass-filter using the Parks-McClellan Algorithm
%
%   Usage: b = pm_filter(order,wpass,wstop)
%
%   Input parameter:
%     order   - order N of filter in original (not upsampled) domain
%     wpass   - last normalised passband frequency [0..1] 
%     wstop   - first normalised stopband frequency [0..1]
%
%   Output parameter:
%     b   - filter coefficients / [(order+1) x 1]
%
%   See also: delayline, thiran_filter

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


%% ===== Computation =====================================================
persistent pmCachedOrder
persistent pmCachedWpass
persistent pmCachedWstop
persistent pmCachedCoefficients

if isempty(pmCachedOrder) || pmCachedOrder ~= order ...
    || isempty(pmCachedWpass) || pmCachedWpass ~= wpass ...
    || isempty(pmCachedWstop) || pmCachedWstop ~= wstop
  
    A = [1 1 0 0];
    f = [0.0 wpass wstop 1.0]; 
    
    pmCachedOrder = order;
    pmCachedWpass = wpass;
    pmCachedWstop = wstop;
    if ~isoctave
        pmCachedCoefficients = firpm(order,f,A).';
    else
        pmCachedCoefficients = remez(order,f,A).';
    end
end
  
b = pmCachedCoefficients;
