function [f,S] = freq_response_wfs_25d(X,xs,src,L,conf)
%FREQ_RESPONSE_WFS_25D simulates the frequency response for 2.5D WFS
%
%   Usage: [f,S] = freq_response_wfs_25d(X,xs,src,L,[conf])
%
%   Input parameters:
%       X           - listener position (m)
%       xs          - position of virtual source (m)
%       src         - source type of the virtual source
%                         'pw' -plane wave
%                         'ps' - point source
%                         'fs' - focused source
%       L           - array length (m)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       f           - corresponding frequency (x) axis
%       S           - simulated frequency response
%
%   FREQ_RESPONSE_WFS_25D(X,xs,src,L,conf) simulates the frequency
%   response of the wave field at the given position X. The wave field is
%   simulated for the given source type (src) using a WFS 2.5 dimensional
%   driving function in the temporal domain.
%
%   References:
%       Spors2009 - Physical and Perceptual Properties of Focused Sources in
%           Wave Field Synthesis (AES127)
%       Spors2010 - Analysis and Improvement of Pre-equalization in
%           2.5-Dimensional Wave Field Synthesis (AES128)
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: wave_field_mono_wfs_25d, wave_field_time_domain_wfs_25d

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
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
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
[X,xs] = position_vector(X,xs);
isargpositivescalar(L);
isargchar(src);

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
% Plotting result
useplot = conf.useplot;


%% ===== Computation ====================================================
% Calculate the wave field in time domain

% Get the position of the loudspeakers
x0 = secondary_source_positions(L,conf);

% Generate frequencies (10^0-10^5)
f = logspace(0,5,500);
% We want only frequencies until f = 20000Hz
idx = find(f>20000,1);
f = f(1:idx);

S = zeros(1,length(f));

% Activity of secondary sources
x0 = secondary_source_selection(x0,xs,src);
nls = size(x0,1);
% Tapering window
win = tapering_window(x0,conf);
% Get the result for all frequencies
for ii = 1:length(f)
    P = 0;
    % Integration over secondary source positions
    for n = 1:nls

        % ================================================================
        % Secondary source model
        % This is the model for the loudspeakers we apply. We use closed cabinet
        % loudspeakers and therefore the 3D Green's function is our model.
        G = point_source(X(1),X(2),x0(n,1:3),f(ii),conf);

        % ================================================================
        % Driving function D(x0,omega)
        D = driving_function_mono_wfs_25d(x0(n,:),xs,src,f(ii),conf);

        % ================================================================
        % Integration
        %              /
        % P(x,omega) = | D(x0,omega) G(x-xs,omega) dx0
        %              /
        %
        % see: Spors2009, Williams1993 p. 36
        %
        P = P + win(n)*D.*G;
    end
    S(ii) = abs(P);
end


% ===== Plotting =========================================================
if(useplot)
    figure; semilogx(f,db(S));
    ylabel('Amplitude (dB)');
    xlabel('Frequency (Hz)');
end
