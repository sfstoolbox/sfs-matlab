function varargout = time_response_line_source(X,xs,conf)
%TIME_RESPONSE_LINE_SOURCE simulates the time response for a line source at
%the given listener position
%
%   Usage: [s,t] = time_response_line_source(X,xs,conf)
%
%   Input parameters:
%       X           - listener position / m
%       xs          - position of line source / m
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       s           - simulated time response
%       t           - corresponding time axis / s
%
%   TIME_RESPONSE_LINE_SOURCE(X,xs,conf) simulates the impulse response of a
%   line source placed at xs at the given virtual microphone position X.
%   The length in samples of the impulse response is given by conf.N.
%   The actual calculation is done via sound_field_imp() and a loop over time t.
%
%   See also: sound_field_imp, freq_response_line_source

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Developers                             *
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
nargmin = 3;
nargmax = 3;
narginchk(nargmin,nargmax);
isargposition(X);
isargxs(xs);
isargstruct(conf);


%% ===== Configuration ==================================================
fs = conf.fs;
N = conf.N;
showprogress = conf.showprogress;
useplot = conf.plot.useplot;


%% ===== Computation ====================================================
% Disable progress bar and plotting for sound_field_imp()
conf.showprogress = false;
conf.plot.useplot = false;
% Get the position of the loudspeaker from line source position.
% NOTE: its directivity [0 -1 0] will be ignored
x0 = [xs 0 -1 0 1];
% Generate time axis
t = (0:N-1)'/fs;
s = zeros(1,length(t));
for ii = 1:length(t)
    if showprogress, progress_bar(ii,length(t)); end
    % Calculate sound field at the listener position
    p = sound_field_imp(X(1),X(2),X(3),x0,'ls',dirac_imp(),t(ii),conf);
    s(ii) = real(p);
end

% Return parameter
if nargout>0, varargout{1}=s; end
if nargout>1, varargout{2}=t; end


%% ===== Plotting ========================================================
if nargout==0 || useplot
    figure;
    figsize(conf.plot.size(1),conf.plot.size(2),conf.plot.size_unit);
    plot(t*1000,s);
    ylabel('amplitude');
    xlabel('time / ms');
end
