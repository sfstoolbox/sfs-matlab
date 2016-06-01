function varargout = time_response_localwfs(X,xs,src,conf)
%TIME_RESPONSE_LOCALWFS simulates the time response for LOCAL WFS at the given
%listener position
%
%   Usage: [s,t] = time_response_localwfs(X,xs,src,conf)
%
%   Input parameters:
%       X           - listener position / m
%       xs          - position of virtual source / m
%       src         - source type of the virtual source
%                         'pw' -plane wave
%                         'ps' - point source
%                         'fs' - focused source
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       s           - simulated time response
%       t           - corresponding time axis / s
%
%   TIME_RESPONSE_LOCALWFS(X,xs,src,conf) simulates the time response of the
%   sound field at the given position X. The sound field is simulated for the
%   given source type (src) using the sound_field_imp_localwfs function.
%
%   See also: sound_field_imp_localwfs, freq_response_localwfs

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
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargposition(X);
isargxs(xs);
isargchar(src);
isargstruct(conf);


%% ===== Configuration ==================================================
% Plotting result
useplot = conf.plot.useplot;
fs = conf.fs;
showprogress = conf.showprogress;
% Disable progress bar for the sound field function
conf.showprogress = 0;
% Disable plotting, otherwise the sound_field_imp fails
conf.plot.useplot = 0;

% Select secondary source type
if strcmp('2D',conf.dimension)
    greens_function = 'ls';
else
    greens_function = 'ps';
end

%% ===== Computation ====================================================
% Get the position of the loudspeakers
x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,xs,src);
x0 = secondary_source_tapering(x0,conf);
% Generate time axis (0-1000 samples)
t = (0:2000)';
s = zeros(1,length(t));
d = driving_function_imp_localwfs(x0,xs,src,conf);
% If desired a cosine shaped pulse instead of the default dirac pulse could be
% used
%d = convolution(d,hann_window(5,5,10));
for ii = 1:length(t)
    if showprogress, progress_bar(ii,length(t)); end
    % Calculate sound field at the listener position
    p = sound_field_imp(X(1),X(2),X(3),x0,greens_function,d,t(ii),conf);
    s(ii) = real(p);
end

% Return parameter
if nargout>0, varargout{1}=s; end
if nargout>1, varargout{2}=t; end


%% ===== Plotting ========================================================
if nargout==0 || useplot
    figure;
    figsize(conf.plot.size(1),conf.plot.size(2),conf.plot.size_unit);
    plot(t/fs*1000,db(abs(s)));
    ylabel('amplitude / dB');
    xlabel('time / ms');
end
