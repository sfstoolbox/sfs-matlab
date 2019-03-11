function [z,p,k] = linkwitz_riley(n,wc,ftype,domain)
%LINKWITZ_RILEY zero-poles-gain representation of the Linkwitz-Riley filter
%
%   Usage: [z,p,k] = linkwitz_riley(n,wc,ftype,[domain])
%
%   Input parameters:
%     n       - order of filter (only even allowed)
%     wc      - cutoff frequency, [0..1] for z-Domain, [0..] rad/s in s-Domain
%     ftype   - filter type {'low','high','all'}
%     domain  - 's' for analog/Laplace-Domain
%               'z' for digital/z-Domain (default)
%                  
%   Output parameters:
%     z       - zeros of filter
%     p       - poles of filter
%     k       - gain of filter
%
%   References:
%       Lipshitz and Vanderkooy (1986), "In-Phase Crossover Network Design",
%       Journal of the Audio Engineering Society, vol. 34, no. 11, pp. 889-894,
%       http://www.aes.org/e-lib/browse.cfm?elib=5237

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2019 SFS Toolbox Developers                             *
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
% https://sfs.readthedocs.io                            sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ===================================
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargpositivescalar(n,wc);
if mod(n,2)
    error('%s: n (%d) is not an even integer.',upper(mfilename),n);
end
isargchar(ftype)
if ~any(strcmp(ftype, {'low', 'high', 'all'}))
    error('%s: ftype (%s) is not a supported filter type.',upper(mfilename), ...
        ftype);
end
if nargin == nargmin
    domain = 'z';
else
    isargchar(domain);
    if ~any(strcmp(domain, {'z', 's'}))
        error('%s: domain (%s) must be either "z" or "s".',upper(mfilename), ...
            domain);
    end
end

%% ===== Computation =====================================================
switch ftype
case  {'low','high'}
    % === lowpass or highpass LR Filter (squared Butterworth Filter) ===
    % See Lipshitz and Vanderkooy (1986), eq. (6) & (12)
    [z,p,k] = butter(n/2,wc,ftype,domain);
    z = [z(:); z(:)];  % octave creates row vectors
    p = [p(:); p(:)];  % octave creates row vectors    
    k = k.^2;
    if strcmp(ftype, 'high')
        k = k.*(-1).^(n/2);
    end
case 'all'
    % === allpass LR Filter (same phase as lowpass and highpass LR Filter) ===
    % See Lipshitz and Vanderkooy (1986), eq. (11)
    [~,p,~] = butter(n/2,wc,'low',domain);
    p = p(:);  % octave creates row vectors
    if strcmp(domain,'z')
        z = 1./conj(p);
        k = prod(p).*(-1).^(n/2);
    else
        z = -p;
        k = (-1).^(n/2);
    end
end

% Average complex conjugates to make them exactly symmetrical.
% This avoids the following Octave bug: https://savannah.gnu.org/bugs/?49996
if ~isempty(z)
    [~, idx] = sort(imag(z));
    z = z(idx);
    z = 0.5*z + 0.5*conj(z(end:-1:1));
end
if ~isempty(p)
    [~, idx] = sort(imag(p));
    p = p(idx);
    p = 0.5*p + 0.5*conj(p(end:-1:1));
end
