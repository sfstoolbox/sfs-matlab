function h = general_least_squares(samples,fractional_delay,passband_edge)
%GENERAL_LEAST_SQUARES fractional delay approximation using the general least
%squares method
%
%   Usage: h = general_least_squares(samples,fractional_delay,passband_edge)
%
%   Input parameters:
%       samples          - filter length (filter order N = samples-1)
%       fractional_delay - fractional delay (0 < x <= 1)
%       passband_edge    - passband edge of approximation (in [0 1])
%
%   Output parameters:
%       h                - filter coefficient vector h(1)...h(samples)
%
%   GENERAL_LEAST_SQUARES(samples,fractional_delay,passband_edge) calculates the
%   filter coefficients needed for the least squares fractional delay method
%   which is used in delayline() if conf.usefracdelay == true and
%   conf.fracdelay_method == 'least_squares'.
%
%   See also: delayline

% This code is based on hgls2() from Timo Laakso

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
% Disabled for time constrains
%nargmin = 3;
%nargmax = 3;
%narginchk(nargmin,nargmax);


%% ===== Computation =====================================================
N = samples-1;                   % filter order
M = N/2;                         % middle value
if (M-round(M))==0
    D = fractional_delay + M;    % integer part closest to middle
else
    D = fractional_delay + M-0.5;
end

cT = zeros(N+1,1);
p1 = cT;
cT(1)=passband_edge;
if round(D)==D
    p1(1) = passband_edge;
else
    % matlab's sinc(x) equals sin(pi*x)./(pi*x)
    p1(1) = passband_edge * sinc(D*passband_edge);
end
for k=1:N  % compute the elements of the Toeplitz matrix (vector)
    k1 = k+1;
    kD = k-D;
    cT(k1) = passband_edge * sinc(k*passband_edge);
    p1(k1) = passband_edge * sinc(kD*passband_edge);
end
P = toeplitz(cT);
h = P\p1;
