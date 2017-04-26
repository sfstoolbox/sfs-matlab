function [zz,pz,kz] = linkwitz_riley(n,wc,ftype)
%LINKWITZ_RILEY computes zero-poles-gain representation in z-domain of
%Linkwitz-Riley filter
%
%   Usage: [zz,pz,kz] = linkwitz_riley(n,wc,ftype)
%
%   Input parameter:
%     n     - order of filter (only even allowed)
%     wc    - normalised cutoff frequency [0..1] 
%     ftype - filter type {'low','high','all'}
%
%   Output parameter:
%     zz    - zeros of filter in z-domain
%     pz    - poles of filter in z-domain
%     kz    - gain of filter

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


%% ===== Checking of input  parameters ========================================
nargmin = 3;
nargmax = 3;
narginchk(nargmin,nargmax);
if mod(n,2)
  error('%s: n (%d) is not an even integer',upper(mfilename),n);
end


%% ===== Configuration ========================================================
switch ftype
case  {'low','high'}
    % === lowpass or highpass LR Filter (squared Butterworth Filter) ===
    [zz,pz,kz] = butter(n/2,wc,ftype);  
    zz = [zz(:); zz(:)];  % octave creates row vectors
    pz = [pz(:); pz(:)];  % octave creates row vectors
    kz = kz.^2;
case 'all'
    % === allpass LR Filter (same phase as lowpass and highpass LR Filter) ===
    [~,pz,~] = butter(n/2,wc,'low');
    pz = pz(:);  % octave creates row vectors
    zz = 1./conj(pz);
    kz = prod(pz);
otherwise
    error('%s: ftype (%s) is not a supported filter type',upper(mfilename), ...
        ftype);
end
