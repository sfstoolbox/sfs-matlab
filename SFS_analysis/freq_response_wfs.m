function varargout = freq_response_wfs(X,xs,src,conf)
%FREQ_RESPONSE_WFS simulates the frequency response for WFS at the given
%listener position
%
%   Usage: [S,f] = freq_response_wfs(X,xs,src,[conf])
%
%   Input parameters:
%       X           - listener position / m
%       xs          - position of virtual source / m
%       src         - source type of the virtual source
%                         'pw' -plane wave
%                         'ps' - point source
%                         'fs' - focused source
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       S           - simulated frequency response
%       f           - corresponding frequency axis / Hz
%
%   FREQ_RESPONSE_WFS(X,xs,src,conf) simulates the frequency response of the
%   wave field at the given position X. The wave field is simulated for the
%   given source type (src) using a monochromatic WFS driving function.
%
%   References:
%       Spors2009 - Physical and Perceptual Properties of Focused Sources in
%           Wave Field Synthesis (AES127)
%       Spors2010 - Analysis and Improvement of Pre-equalization in
%           2.5-Dimensional Wave Field Synthesis (AES128)
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: wave_field_mono_wfs, wave_field_imp_wfs

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox      sfstoolbox@gmail.com  *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargposition(X);
isargxs(xs);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
end
isargstruct(conf);


%% ===== Configuration ==================================================
% Plotting result
useplot = conf.plot.useplot;


%% ===== Computation ====================================================
% Get the position of the loudspeakers
x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,xs,src);
x0 = secondary_source_tapering(x0,conf);
% Generate frequencies (10^0-10^5)
f = logspace(0,5,500);
% We want only frequencies until f = 20000Hz
idx = find(f>20000,1);
f = f(1:idx);
S = zeros(1,length(f));
nls = size(x0,1);
% Get the result for all frequencies
for ii = 1:length(f)
    % calculate wave field at the listener position
    P = wave_field_mono_wfs(X(1),X(2),X(3),xs,src,f(ii),conf);
    S(ii) = abs(P);
end

% return parameter
if nargout>0 varargout{1}=S; end
if nargout>1 varargout{2}=f; end

% ===== Plotting =========================================================
if nargout==0 || useplot
    figure; semilogx(f,db(S));
    ylabel('Amplitude (dB)');
    xlabel('Frequency (Hz)');
end
